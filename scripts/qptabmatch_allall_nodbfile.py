#!/usr/bin/env python
###############################################################################
#
# qptabmatch_allall_nodbfile.py - 
#                        run the QP tableau matching with all tableauxdistmatrix
#                        files in a directory against each other
#
# File:    qptabmatch_allall_nodbfile.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Run QP tableau match on all tableaux in a directory against each other
# Each file in the directory
# has a .tableaudistmatrix suffix and contains the header line with
# identifeir and dimension, then tabelau and SSE distance matrix
# (both lower triangle fortran format
# for tsrchrd_sparse etc.)
#
# This is _allall_nodbfile since for Fischer data set the db is actually just
# all the tableaux in directory, so it is doing
# duplicate comparisons so for n (=68) tableux it does
# n*n (=4624) comparisons.
#
# Note _allall uses a dbfile so n comparisons are done each run of tsrchd,
# so much more efficient (less overhead) than doing it this way.
# This _allall_nodbfile version does only one comparison on each run,
# the way MSVNS4MaxCMO (Pelta et al 2008) does, to make timings comparable.
# 
#
# Usage:
#     qptabmatch_allall_nodbfile.py query_directory results_directory
#
# query_directory is the directory containing .tableaudistmatrix files,
# as built with build_fischer_db.sh for example.
#
# results_dirctory is a directory to write the output to.
# Each query (.tableauxdistmatrix file) results in one file created
# with .out suffix
# in the output directory, containing the results from that query
# against all the other .tableauxdistmatrix files (and itself).
# Each file is created by parsing  tsrchd_sparse output:
# each line is query identifier and score (whitespace delimited).
# WARNING: these files overwritten if they exist.
# results_directory is created if it does not exist.
#
# Environment variables:
#
#   PATH must contain the location of tsrchd_sparse.
#
#
# $Id: qptabmatch_allall_nodbfile.py 1970 2008-10-10 01:12:47Z astivala $
# 
###############################################################################

import sys,os,glob

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <query_directory> <results_directory>\n")
    sys.exit(1)

    
def main():
    """
    main for qptabmatch_allall_nodbfile.py
    """
    if len(sys.argv) != 3:
        usage(os.path.basename(sys.argv[0]))

    query_directory = sys.argv[1]
    results_directory = sys.argv[2]

    if not os.path.exists(results_directory):
        os.mkdir(results_directory)
    elif not os.path.isdir(results_directory):
        sys.stderr.write('%s is not a directory\n' % results_directory)
        sys.exit(1)

    input_list = glob.glob(os.path.join(query_directory, '*.tableaudistmatrix'))
    i = 0
    while i < len(input_list):
        qfile = input_list[i]
        qid = open(qfile).readline()[:8].lstrip().rstrip()
        if qid == '':
            i += 1
            continue  # empty file, happens for 41010
        outfile = os.path.join(results_directory, 
                        os.path.splitext(os.path.basename(qfile))[0] + '.out' )
        outfh = open(outfile, 'w')
        j = 0
        while j < len(input_list):
            tgtfile = input_list[j]
            if os.path.getsize(tgtfile) == 0:
                j += 1
                continue # empty file, happens for 41010
            (tsrchd_in, tsrchd_out) = os.popen2('tsrchd_sparse')
            tsrchd_in.write(tgtfile + '\n')   # name of db file
            tsrchd_in.write('T T F\n') # LTYPE LORDER LSOLN
            tsrchd_in.write(open(qfile).read()) # tableau+distmatrix of qfile
            tsrchd_in.close()
            for line in tsrchd_out:
                if line[0] == '#':
                    continue
                sline = line.split()
                if sline[0] == '':
                    continue
                tgtid = sline[0]
                score = float(sline[1])
                outfh.write('%8s %12.4f\n' % (tgtid,score))
            j += 1
        outfh.close()
        i += 1

if __name__ == "__main__":
    main()

