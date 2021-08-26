#!/usr/bin/env python
###############################################################################
#
# superimposessemap.py - Superimpose structures according to SSE mapping
#
# File:    superimposessemap.py
# Author:  Alex Stivala
# Created: August 2008
#
#
# Supermipose in 3D the residues in corresponding SSEs by orthogonal
# transformations (using SVD) using the Bio.PDB.Superimposer module.
#
# $Id: superimposessemap.py 1821 2008-08-18 00:54:56Z astivala $
# 
###############################################################################

"""

Using the SSE mapping from soln2ssemap.py, which shows pairs of SSE
sequential (from 1) numbers that correspond to each other, use orthogonal
transormation to superimpose the residues in corresponding SSEs,
calculating RMSD and producing superimposition in a PDB file for visualization.

Requires the ptsecstruct.py module to get secondary structures using
DSSP (or STRIDE) (add directory contianing to to PYTHONPATH).

Note that these must be the same definintions used
to produce the mapping, i.e. that the tableaux database and query
were built with, otherwise it won't realy make sense.

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
from time import strftime,localtime

import Bio.PDB

import ptsecstruct
from ptutils import biopdbresid_to_pdbresseq,get_int_icode

from parsessemap import parse_ssemap,SearchMap,QuerySSEMap
from pathdefs import ASTRAL_ROOT



#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def get_structure(scopsid, thepdbfile=None):
    """
    Get Bio.PDB parsed structure for specified identifier or PDB file.

    Parameters:
       scopsid - SCOP identifier to get SSEs for; used to locate file
                 under ASTRAL SCOP hierarchy.
       thepdbfile - (default None) if not None, PDB file to get SSEs for,
                 overriding scopsid.
    Return value:
       Bio.PDB parsed structure.
    """
    if thepdbfile:
        pdbfile = thepdbfile
    else:
        pdbfilename = os.path.join(scopsid[2:4].lower(),
                                   scopsid.lower() + '.ent')
        pdbfile = os.path.join(ASTRAL_ROOT, pdbfilename)

    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(scopsid, pdbfile)
    return structure


def get_sse_nodes(scopsid, thepdbfile=None):
    """
    Get SSE definitions in form of PTNode objects
    from the supplied SCOP sid using
    DSSP. Uses the ptsecstruct.py module, note comments at top of this
    module also regarding ensuring the same definitions are used here
    as for the actual search.

    Parameters:
       scopsid - SCOP identifier to get SSEs for; used to locate file
                 under ASTRAL SCOP hierarchy.
       thepdbfile - (default None) if not None, PDB file to get SSEs for,
                 overriding scopsid.
    Return value:
       list of PTNode objects represneting the SSEs.
    """
    if thepdbfile:
        pdbfile = thepdbfile
    else:
        pdbfilename = os.path.join(scopsid[2:4].lower(),
                                   scopsid.lower() + '.ent')
        pdbfile = os.path.join(ASTRAL_ROOT, pdbfilename)

    secstruct = ptsecstruct.read_secstruct_from_dssp(pdbfile)
    return secstruct.get_sse_tuple_list()
    

def get_residue_list(model):
    """
    Get list of Bio.PDB.Residue objects in supplied Bio.PDB.Model
    Parmeters:
       model - Bio.PDB.Model object
    Return value:
       List of Bio.PDB.Residue objects in the model
       
    """
    residue_list = []
    for chain in model:
        # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
        residue_list += [ residue for residue in chain.get_unpacked_list()
                          if Bio.PDB.is_aa(residue) ]
    return residue_list


def build_resid_dict(residue_list):
    """
    Build dictionary mapping (chainid, pdb_resid) to index in residue_list
    for all residues, not just those in this domain.

    Parameters:
        residue_list - list of Bio.PDB.Residue objects
    Return value:
         dict of { {chainid,pdb_resseq) : seqindx }
         where chainid and pdb_resseq make up
         the PDB residue identifier, the pdb_resseq
         being string resnum+icode if any e.g.
         '60' or '60A', seqindx is the indiex
         into sequential list of all residues
         residue_list.
    """
    pdb_resid_dict = {}
    seq_indx = 0
    while seq_indx < len(residue_list):
        residue = residue_list[seq_indx]
        pdb_resid_dict[( ptsecstruct.pdb_chainid_to_stride_chainid(
                         residue.get_full_id()[2]), 
                         biopdbresid_to_pdbresseq(residue.get_id()) )] = seq_indx
        seq_indx += 1
    return pdb_resid_dict

    
def get_matched_residues(matched_sses, query_struct, db_struct):
    """
    Given the list of correpsonding SSEs in the two structures, return
    list of corresponding Bio.PDB.Residue objects.

    Parameters:
       matched_sses - list of (A,B) tuples where A and B are
                      tuples (chain, start_resi, end_resi, type) in
                      query_struct and db_struct respectively.
       query_struct - Bio.PDB.Structure
       db_struct -    Bio.PDB.Structure
    Return value:
       tuple (match_query_residues, match_db_residues) of equal length lists of
       corresponding Bio.PDB.Residue objects in query and db structs resp.
    """
    query_model = query_struct[0] # always using model 0 (TODO)
    db_model = db_struct[0]       # always using model 0 (TODO)

    query_residue_list = get_residue_list(query_model)
    query_resid_dict = build_resid_dict(query_residue_list)
    db_residue_list = get_residue_list(db_model)
    db_resid_dict = build_resid_dict(db_residue_list)

    match_query_residues = []
    match_db_residues = []
    for ((qchain, qstart_resi, qend_resi, qtype),
         (dchain, dstart_resi, dend_resi, dtype)) in matched_sses:
        try:
            start_indx = query_resid_dict[(qchain, qstart_resi)]
        except KeyError:
            # May be HETATM
            while not query_resid_dict.has_key((qchain, qstart_resi)):
                qstart_resi = str(get_int_icode(qstart_resi)[0] + 1)
            start_indx = query_resid_dict[(qchain, qstart_resi)]                
        try:
            end_indx = query_resid_dict[(qchain, qend_resi)]
        except KeyError:
            # May be HETATM
            while not query_resid_dict.has_key((qchain, qend_resi)):
                qend_resi = str(get_int_icode(qend_resi)[0] - 1)
            end_indx = query_resid_dict[(qchain, qend_resi)]
        query_residues = query_residue_list[start_indx : end_indx + 1]
        try:
            start_indx = db_resid_dict[(dchain, dstart_resi)]
        except KeyError:
            # May be HETATM
            while not db_resid_dict.has_key((dchain, dstart_resi)):
                dstart_resi = str(get_int_icode(dstart_resi)[0] + 1)
            start_indx = db_resid_dict[(dchain, dstart_resi)]
        try:
            end_indx = db_resid_dict[(dchain, dend_resi)]
        except KeyError:
            # May be HETATM
            while not db_resid_dict.has_key((dchain, dend_resi)):
                dend_resi = str(get_int_icode(dend_resi)[0] - 1)
            end_indx = db_resid_dict[(dchain, dend_resi)]                
        db_residues = db_residue_list[start_indx : end_indx + 1]

#         # if the SSEs are of unequal length, just truncate the longer
#         # FIXME: should do something better here, e.g. use residues
#         # in middle of SSEs since definitions at ends probably less certain

#         if len(db_residues) > len(query_residues):
#             db_residues = db_residues[:len(query_residues)]
#         elif len(query_residues) > len(db_residues):
#             query_residues = query_residues[:len(db_residues)]

#         match_query_residues += query_residues
#         match_db_residues += db_residues


#         # use the first and last residues in each SSE
#         # FIXME: should really use projected enpoints on vector
#         # to represent the vector actually used to construct tableau
#         # as per fit_axis in ptnode.py
#         match_query_residues += [query_residues[0], query_residues[-1]]
#         match_db_residues += [db_residues[0], db_residues[-1]]


        # another dodgy way: just the 'most cetnral' residue (FIXME)
        match_query_residues.append(query_residues[len(query_residues)/2])
        match_db_residues.append(db_residues[len(db_residues)/2])
    
                                
    assert(len(match_query_residues) == len(match_db_residues))
    return (match_query_residues, match_db_residues)


    
#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-d domainid] [-u query_pdbfile] [-b db_pdbfile] [-o outputdir] \n")
    sys.stderr.write(
        "-d domainid: use this structure, if more than one in input\n"
        "-u query_pdbfile:  filename of query PDB file. If not specified then\n"
        "               identifier is used to find in ASTRAL SCOP hierarchy.\n"
        "-b db_pdbfile:  filename of database PDB file. If not specfied then\n"
        "              identifier is used to find in ASTRAL SCOP hierarchy.\n"
        "           Only valid if there is only one domain (either becuase -d is\n"
        "           specified or there is only one in the input).\n"
        "-o outputdir: directory to write PDB of superimposed structures in.\n"
        )
    sys.exit(1)

    
def main():
    """
    main for superimposessemap.py

    Usage: superimposessemap.py [-d domainid] [-s] [-u query_pdbfile] [-b db_pdbfile] [-o outputdir]

    -d domainid: only output for this domain, not all
    -u query_pdbfile:  filename of query PDB file. If not specified then
                       identifier is used to find in ASTRAL SCOP hierarchy.
    -b db_pdbfile:  filename of database PDB file. If not specfied then
                      identifier is used to find in ASTRAL SCOP hierarchy.
                   Only valid if there is only one domain (either becuase -d is
                   specified or there is only one in the input).
    -o outputdir: directory to write PDB files of superimposed structures in.

                 
    Input is on stdin, the output of soln2ssemap.py,
    identifier and score (as per input), then
    for each matching a line containing
    i and j separated by a space,
    one per line (with blank line before next id) e.g.:
    
    d1wiua_     -23.0000
    1 1
    3 2
    8 4
    9 5
    11 6
    14 9


    The first SSE number on each line is in the query structure
    (specified in header information), the second
    is in the db structure (d1wiua_ in example).

    Output is RMSD value on stdout, and PDB file(s) in specified directory if -o
    specfied.
    stdout output format is one result per line, fields whitespace delimited:
    
    identifier score num_sses_matched num_aligned_points rmsd

    e.g.

    d1t10a_ -40.9999 8 16 16.93

    num_aligned_points is number of points used in the superposition,
    RMSD is the RMS deviation of those points (in Angstroms).
    
    """
    global verbose
    verbose = False

    dbdomid = None
    query_pdbfile = None
    db_pdbfile = None
    outputdir = None
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "d:u:b:o:")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-d": # domain id specified, only get this one
            dbdomid = arg
        elif opt == "-u":  # query PDB filename
            query_pdbfile = arg
        elif opt == "-b":  # db PDB filename
            db_pdbfile = arg
        elif opt == "-o":  # output directory 
            outputdir = arg
        else:
            usage(os.path.basename(sys.argv[0]))
        
    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))

    search_maps = parse_ssemap(sys.stdin)

    if (db_pdbfile and not dbdomid and len(search_maps.query_ssemap_list) > 1):
        sys.stderr.write("ERROR: -b specified without -d and more than one "
                         "structure on input\n")
        sys.exit(1)

    query_sse_nodes = get_sse_nodes(search_maps.queryid, query_pdbfile)
    query_structure = get_structure(search_maps.queryid, query_pdbfile)
    for query_ssemap in search_maps.query_ssemap_list:
        if ((not dbdomid) or (query_ssemap.domid == dbdomid)):
            db_sse_nodes = get_sse_nodes(query_ssemap.domid, db_pdbfile)
            db_structure = get_structure(query_ssemap.domid, db_pdbfile)
            sse_map = query_ssemap.sse_map
            if len(sse_map) == 0:
                sys.stderr.write('no SSEs matched for ' + query_ssemap.domid +
                                 ': skipping\n')
                continue
            
            matched_sse_nodes = [(query_sse_nodes[i-1],db_sse_nodes[j-1]) for (i,j) in sse_map]
            matched_residues = get_matched_residues(matched_sse_nodes,
                                                    query_structure,
                                                    db_structure)

            # get Carbon alpha atoms for matched residues
            query_atoms = [residue['CA'] for residue in matched_residues[0]]
            db_atoms = [residue['CA'] for residue in matched_residues[1]]

            # get orthogonal transformation to superimpose query and db atoms
            superimposer = Bio.PDB.Superimposer()
            superimposer.set_atoms(query_atoms, db_atoms)

            # get the RMSD for the atoms used to calculate transformation
            rmsd = superimposer.rms

            sys.stdout.write('%s %8.4f %4d %4d %6.2f\n' %
                             (query_ssemap.domid,query_ssemap.score,
                              len(sse_map),
                              len(matched_residues[0]), rmsd))

            if outputdir:
                if not os.path.isdir(outputdir):
                    sys.stderr.write("'" + outputdir + "' is not an existing "
                                     "directory, no output written\n")
                else:
                    # apply the transformation to all db atoms
                    superimposer.apply(db_structure.get_atoms())
                    
                    # save aligned structure as PDB file
                    io = Bio.PDB.PDBIO()
                    io.set_structure(db_structure)
                    outpdbfilename = search_maps.queryid.lstrip().rstrip() + \
                                     '_' + \
                                     query_ssemap.domid.lstrip().rstrip() + \
                                     '.pdb'
                    outpdbfh = open(os.path.join(outputdir,outpdbfilename), 'w')
                    outpdbfh.write('REMARK generated by ' +
                                   os.path.basename(sys.argv[0]) + '\n')
                    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
                    outpdbfh.write('REMARK on ' + timestamp + '\n')
                    outpdbfh.write('REMARK \n')
                    outpdbfh.write('REMARK ' + query_ssemap.domid +
                                   ' superimposed on ' + search_maps.queryid +
                                   '\n')
                    outpdbfh.write('REMARK SCORE = %8.4f\n' % query_ssemap.score)
                    outpdbfh.write('REMARK NSSES = %4d\n' % len(sse_map))
                    outpdbfh.write('REMARK NRES = %4d\n' % len(matched_residues[0]))
                    outpdbfh.write('REMARK RMSD = %6.2f\n' % rmsd)
                    outpdbfh.write('REMARK \n')
                    outpdbfh.write('REMARK from:\n')
                    for cline in search_maps.comment_lines:
                        outline = cline[:65]
                        outpdbfh.write('REMARK ' + outline)
                        if outline[-1] != '\n':
                            outpdbfh.write('\n')
                    io.save(outpdbfh)
                    outpdbfh.close()
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
