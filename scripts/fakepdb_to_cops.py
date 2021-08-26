#!/usr/bin/env python
#
# File:    fakepdb_to_cops.py
# Author:  Alex Stivala
# Created: April 2010
#
# $Id: fakepdb_to_cops.py 3909 2010-07-11 04:45:25Z alexs $
"""
 fakepdb_to_cops.py - Convert fake PDB identifers back to COPS identifiers

 Usage: fakepdb_to_cops.py [fakepdbids_filename] < DaliteLite-2col-output

 The input file is 2 column from dalilitout2col.py
 on stdin.
 Output is to stdout.

  If the fakepdbids_filename is supplied, the translation table is read
  from it, toerhwise the COPS file is used.

 Note that DaliLite ONLY allows 4 char PDB codes, with chain appended, so
 SCOP or COPS type codes will not work in cases where there are more than
 one with the same pdb code and chain e.g. c1d3uB1 and c1d3uB2
 This is highly invonvenient, we have to get around it by mapping ALL
 structures to 'fake' PDB codes with the cops_to_fakepdb.py and back with the
 fakepdb_to_cops.py scripts.

"""

import sys,os

from tsevalutils import iter_searchresult


COPS_FAKEPDBIDS_FILE = "/home/alexs/phd/qptabsearch/data/COPS/cops.fakepdbids"



def parse_fakepdbids_file(fname):
    """
    Parse the fake pdb id file to build dictionary mapping fake pdb id
    to COPS Id

    Parameters:
       fname - name of file to parse, first col is COPS id, 2nd is fake PDB id

    Return value:
       tuple ( fake2cops, cops2fake ) where
        fake2cops is dict { fakepdb : copsid } mapping fake PDB id to COPS id
        cops2fake is dict { copsid  : fakepdb } mapping COPS to fake PDB id
    """
    fake2cops_dict = {}
    cops2fake_dict = {}
    for line in open(fname):
        if line[0] == '#':
            continue
        sline = line.split()
        if len(sline) != 2:
            sys.stderr.write('bad line: %s\n' % line)
        fake2cops_dict[sline[1]] = sline[0]
        cops2fake_dict[sline[0]] = sline[1]
    return ( fake2cops_dict, cops2fake_dict )


def usage(progname):
    sys.stderr.write("Usage: " + progname + " < DaliLite2ColOut\n")
    sys.exit(1)


def main():
    """
    main for fakepdb_to_cops.py - see usage message at file header
    """
    if len(sys.argv) == 2:
        fakeids_filename = sys.argv[1]
    elif len(sys.argv) == 1:
        fakeids_filename = COPS_FAKEPDBIDS_FILE
    else:
        usage(os.path.basename(sys.argv[0]))

    FAKE_TO_COPS_DICT = parse_fakepdbids_file(fakeids_filename)[0]

    for (score, dbid) in iter_searchresult(sys.stdin,multiquery=False):
        sys.stdout.write("%s    %f\n" % (FAKE_TO_COPS_DICT[dbid[:4]], score))

            
if __name__ == "__main__":
    main()
