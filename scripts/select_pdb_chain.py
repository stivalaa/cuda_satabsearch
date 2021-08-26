#!/usr/bin/env python
###############################################################################
#
# select_pdb_chain.py - select just one chain from a PDB file
#
#
# File:    select_pdb_chain.py
# Author:  Alex Stivala
# Created: April 2010
#
#
# $Id: select_pdb_chain.py 3572 2010-04-20 04:07:48Z alexs $
#
###############################################################################
"""
This script parses a PDB file and writes out another PDB file with
only the specified chain in it, using BioPython.

See usage in main() documentation.
"""

import os,sys
import gzip
from Bio.PDB import *

class ChainSelect(Select):
    """
    The ChainSelect class inherits from the PDBIO.Select class
    and overrides function to select only certain residues for writing
    ATOM records in the chain we are interested in.

    See the Bio.PDB documentation by Thomas Hamelryck:
      biopython-1.43/Doc/biopdb_faq.pdf
    """
    def __init__(self, chainid):
        """
        Constructor for the ChainSelect class, sets the chainid member
        used to accept only residues in that chain.
        Parameters:
           chainid - chain id to select
        """
        self.chainid = chainid.upper()

    def __repr__(self):
        """
        Overrides the base __repr__ to write out the domain we have
        """
        return "<ChainSelect: " + self.chainid + ">"
    

    def accept_chain(self, chain):
        """
        Overrides the base accept_chain() to select only the chain we want

        Parameters:
            chain - Bio.PDB Chain object for the chain in question

        Return value:
            1 to accept the chain (its chainid is the one we want)
            0 otherwise (do not select the chain)
        """
        if chain.get_id().upper() == self.chainid:
            return 1
        else:
            return 0


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " pdbfile chainid\n")
    sys.exit(1)


def main():
    """
    main for select_pdb_chain.py

    Usage:
        select_pdb_chain.py pdbfile chainid

    Parses the specified PDB in pdbfile and writes out a version with
    only the specified chain by chainid in it.
    
    Filenames may be either in the 1QLP.pdb format or the pdbq1lp.ent format.
    Compressed pdb files are supported (gzip) (e.g. pdb1qlp.ent.gz).

    The output filename is the input with the chain appened after an
    undserscore e.g. 1QLP_A.pdb.
    """

    if len(sys.argv) != 3:
        usage(os.path.basename(sys.argv[0]))


    pdb_filename = sys.argv[1]
    chainid = sys.argv[2]

    if len(chainid) != 1:
        usage(os.path.basename(sys.argv[0]))
        
    pdb_file_basename = os.path.basename(pdb_filename)
    (name,extension) = os.path.splitext(pdb_file_basename)

    if extension == '.gz':
        pdb_fh = gzip.open(pdb_filename)
    else:
        pdb_fh = open(pdb_filename)

    if len(name) == 4:
        pdbid = name
    elif len(name) >= 7 and name[:3] == 'pdb':
        pdbid = name[3:7]

    outfilename = pdbid + '_' + chainid + '.pdb'

    parser = PDBParser()
    structure = parser.get_structure(pdbid, pdb_fh)
    pdb_fh.close()

    io = PDBIO()
    io.set_structure(structure)
    io.save(outfilename, ChainSelect(chainid))


if __name__ == "__main__":
    main()
