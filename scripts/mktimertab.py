#!/usr/bin/env python
#
# File:    mktimertab.py
# Author:  Alex Stivala
# Created: November 2009
#
#
# mktimertab.sh - make table for R read.table(,header=TRUE) with CPU times
#                 from tsrchd_* -t output
#  
#
# Usage: mktimertab.sh  tsrdcd-t_outputfilename
#
# The input file is the output from tsrchd_* -t (3 columns:
# queryid, score, cputime), with QUERY ID and DBFILE In comments)
# from which we extract the scopsid as the first component (e.g. d1ubia_)
# and use  ../data/queryid.input and the database file 
#
# Output is to stdout, with each row containing
#
#  queryid dbid querysses dbsses score cputime
#
#
# (Note that queryid and querysses will actually be the same for every row)
#
# $Id: mktimertab.py 2898 2009-11-02 04:14:06Z astivala $
#


import os,sys

# location of .input files (used for getting number of SSEs and db filename)
INPUTDIR = os.getenv("HOME") + "/phd/qptabsearch/data"

if len(sys.argv) != 2:
    sys.stderr.write("Usage: " + sys.argv[0] + " tsrchd-t_filename\n")
    sys.exit(1)

tsrchdfile = sys.argv[1]


sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
sys.stdout.write('# results from:\n')
firsttime = True
for line in open(tsrchdfile):
    if line[:12] == "# QUERY ID =":
        sys.stdout.write("# " + line)
        queryid = line.split("=")[1].lstrip().rstrip().lower()
        inputfile = os.path.join(INPUTDIR, queryid + ".input")
        for infline in open(inputfile):
            if infline[:len(queryid)].lower() == queryid:
                querysses = infline.split()[1]
                break
    elif line[:10] == "# DBFILE =":
        sys.stdout.write("# " + line)
        dbfile = line.split("=")[1].lstrip().rstrip()
    elif line[0] == "#" or len(line) == 0:
        sys.stdout.write("# " + line)
    else:
        if firsttime:
            sys.stdout.write("queryid dbid querysses dbsses score cputime\n")
            firsttime = False
            # build dict of number of SSEs in each database structure
            dbnumsse_dict = {} # dict of {scopid : numsses}
            for dbline in open(dbfile):
                if dbline[0] == "d":
                    splitdbline = dbline.split()
                    dbnumsse_dict[splitdbline[0]] = splitdbline[1]
        splitline = line.split()
        dbid = splitline[0]
        score = splitline[1]
        cputime = splitline[2]
        dbsses = dbnumsse_dict[dbid]
        sys.stdout.write("%s %s %s %s %s %s\n"
                         % (queryid, dbid, querysses, dbsses, score, cputime))
