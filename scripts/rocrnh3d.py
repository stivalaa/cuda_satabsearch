#!/usr/bin/env python
#
#
# rocrnh3d.py - Output score+label table for results on NH3D data set
#
# File:    rocnh3d.py
# Author:  Alex Stivala
# Created: September 2008
#
# Evaluate QP tableau search or MSVNS4MaxCMO for all against all matching
# for the Nh3D data set (Thiruv et al 2005 BMC Struct Biol 5:12)
# as per Pelta et all 2008 BMC Bioinformatics 9:161
#
# $Id: rocrnh3d.py 2079 2009-03-03 07:43:11Z astivala $
# 

"""
Write scores from structure search method and actual class labels
(0/1 for same/different fold as query). Using  CATH architecture as gold 
standard.
Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os,glob
import getopt

from cathmap import CATHMAP
from tsevalutils import compute_auc,parse_searchresult

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# list of different CATH architectures 
ARCH_LIST= ["1.10", "1.20", "2.10", "2.170", "2.30", "2.40", "2.60", "2.70", "3.10", "3.20", "3.30", "3.40", "3.60", "3.90", "4.10"]

# List of query CATH identifiers, from the Additional File 1 spreadsheet
# for Pelta et al 2008
QUERY_LIST=["1.10.1040", "1.10.1320", "1.10.533", "1.10.645", "1.20.1280", "1.20.210", "1.20.5", "1.20.840", "2.10.25", "2.10.260", "2.10.270", "2.10.90", "2.170.16", "2.170.230", "2.170.290", "2.170.40", "2.30.110", "2.30.18", "2.30.230", "2.30.29", "2.30.40", "2.40.155", "2.40.160", "2.40.180", "2.40.340", "2.40.50", "2.60.130", "2.60.260", "2.60.420", "2.60.90", "2.70.100", "2.70.180", "2.70.220", "2.70.98", "3.10.105", "3.10.170", "3.10.270", "3.10.330", "3.10.400", "3.20.120", "3.20.140", "3.20.19", "3.20.70", "3.20.90", "3.30.1530", "3.30.1690", "3.30.240", "3.30.559", "3.30.560", "3.30.60", "3.30.990", "3.40.1210", "3.40.1380", "3.40.225", "3.40.720", "3.60.100", "3.60.120", "3.60.20", "3.60.40", "3.60.90", "3.90.1280", "3.90.1300", "3.90.1350", "3.90.1580", "3.90.510", "3.90.850", "4.10.1080", "4.10.1090", "4.10.220", "4.10.260", "4.10.480", "4.10.540", "4.10.790"]

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-cnv] "
                     "  <outdir>\n")
    sys.stderr.write('  -c evaulate at Class not Architecture level\n')
    sys.stderr.write('  -n negate scores\n')
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.exit(1)

    
def main():
    """
    main for rocnh3d.py

    Usage: rocnh3d.py [-cnv]   <outdir>

    
    -c evaluate at CATH Class not Architecture level
    -v turns on debug output to stderr
    -n negate scores

    <outdir> is the directory containing output files as generated
    by tsrchd_sparse (qptabmatch_allall.py) or msvns4maxcmo_allall.py 

    The table of scores and class labels is printed to stdout.
    """
    global verbose
    verbose = False
    negateflag = False
    class_level = False
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "cvn?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        elif opt == "-c": #class not arch
            class_level = True
        elif opt == "-n": # negate scores
            negateflag = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    outdir = args[0]

    # build dict of search results, one for each query
    # and corresponding dict of gold standard results
    searchresult_dict = {}
    goldstd_dict = {}
    for result_file in glob.glob(os.path.join(outdir, '*.out')):
        query_id = os.path.splitext(os.path.basename(result_file))[0].lower()
        result_fh = open(result_file)
        (searchresult,commentlist) = parse_searchresult(result_fh, negateflag)
        result_fh.close()
        searchresult_dict[query_id] = searchresult

        # get the gold standard as the list of 'compreseed' CATH ids
        # that have same architecture as query
        query_id_split = query_id.split('.')
        query_class = query_id_split[0]
        query_arch = query_id_split[1]
        goldstd_ids = []
        for (compressed_id, cathid) in CATHMAP.iteritems():
            cathid_split = cathid.split('.')
            cathid_class = cathid_split[0]
            cathid_arch = cathid_split[1]
            if cathid_class == query_class:
                # only check Architecture match if not evaluating at Class level
                if class_level or cathid_arch == query_arch:
                    goldstd_ids.append(compressed_id)

        goldstd_dict[query_id] = goldstd_ids


    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us

    sys.stdout.write('score    label\n')
    sys.stdout.write('#-------------\n')
    for query_id in searchresult_dict.iterkeys():
        searchresult = searchresult_dict[query_id]
        goldstd_pos_dict = dict([(domainid,True) for domainid in 
                                 goldstd_dict[query_id]])
        for (score, domainid) in searchresult:
            if goldstd_pos_dict.has_key(domainid):
                label = 1
            else:
                label = 0
            sys.stdout.write('%20.8f  %d\n' % (score,label))
    
    
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()

