#!/usr/bin/env python
###############################################################################
#
# tsevalfischer.py - Evaluate tableaux search or MAX-CMO 
#                    against Fischer data set
#
# File:    tsevalfischer.py
# Author:  Alex Stivala
# Created: September 2008
#
# Evaluate QP tableau search or MSVNS for MAX-CMO (Pelta et al 2008)
# for all against all matching
# for the Fischer data set (Fischer et al 1996 Pac. Symp. Biocomput. 300-318))
# as per Pelta et all 2008 BMC Bioinformatics 9:161
#
# $Id: tsevalfischer.py 1956 2008-10-07 00:28:16Z astivala $
# 
###############################################################################

# OBSOLETE - use rocrfischer.py and rocauc.r now

"""
Evaluate false negatives using Fischer Table II at each cutoff (rank/score),
domains that are not included above the cuttoff but are in same fold 
(for fold level evaluation) or class (for class level evaluation)
in Table II of Fischer 1996 are false negatives.

Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt

from tsevalutils import compute_auc,parse_searchresult,eval_fn

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# The 68 probe sequences from Fischer 1996 Table II
# Note several PDB ids obsoleted, so change to the replacments

# map id to name of fold
FISCHER_ID_FOLD_DICT = {
    '1dxt_b' : 'globin-like',
    '1cpc_l' : 'globin-like',
    '1c2r_a' : 'cytochrome',
    '2mta_c' : 'cytochrome',
    '1bbh_a' : 'helical bundle',
    '1bge_b' : 'helical bundle',
    '1rcb'   : 'helical bundle',
    '1aep'   : 'helical bundle',
    '1osa'   : 'ef-hand',
    '2sas'   : 'ef-hand',
    '1hom'   : 'other alpha',
    '1lga_a' : 'other alpha',
    '2hpd_a' : 'other alpha',
    '1chr_a' : 'tim barrel',
    '2mnr'   : 'tim barrel',
    '3rub_l' : 'tim barrel',
    '1crl'   : 'hydrolase',
    '1tah_a' : 'hydrolase',
    '1aba'   : 'thieredoxin',
    '1dsb_a' : 'thieredoxin',
    '1gpl_a' : 'thieredoxin',
    '1atn_a' : 'ribonuclease',
    '1hrh_a' : 'ribonuclease',
    '3chy'   : 'open sheet',
    '2ak3_a' : 'open sheet',
    '1gky'   : 'open sheet',
    '2cmd'   : 'open sheet',
    '1eaf'   : 'open sheet',
    '2gbp'   : 'open sheet',
    '1mio_c' : 'open sheet',
    '2pia'   : 'open sheet',
    '1gal'   : 'open sheet',
    '1npx'   : 'open sheet',
    '2hhm_a' : 'mixed',
    '1hip'   : 'small',
    '1isu_a' : 'small',
    '1fc1_a' : 'ig',
    '2fbj_l' : 'ig',
    '1cid'   : 'ig-like',
    '1pfc'   : 'ig-like',
    '1ten'   : 'ig-like',
    '1tlk'   : 'ig-like',
    '3cd4'   : 'ig-like',
    '3hla_b' : 'ig-like',
    '1aaj'   : 'copredoxin',
    '2afn_a' : 'copredoxin',
    '2aza_a' : 'copredoxin',
    '4sbv_a' : 'virus',
    '1bbt_1' : 'virus',
    '1sac_a' : 'lectin-like',
    '1lts_d' : 'ob-fold',
    '1tie'   : 'trefoil',
    '8i1b'   : 'trefoil',
    '1arb'   : 'trypsin',
    '2sga'   : 'trypsin',
    '2snv'   : 'trypsin',
    '1mdc'   : 'lipocalin',
    '1mup'   : 'lipocalin',
    '2sim'   : 'propeller',
    '1cau_b' : 'other beta',
    '2omf'   : 'other beta',
    '1fxi_a' : 'ub fold',
    '1cew'   : 'cystatin',
    '1stf_i' : 'cystatin',
    '2pna'   : 'sh2',
    '2sar_a' : 'other alpha+beta',
    '1onc'   : 'other alpha+beta',
    '5fd1'   : 'other alpha+beta'
}

# map name of fold to list of ids
FISCHER_FOLD_IDLIST_DICT = {
    'globin-like'    : ['1dxt_b','1cpc_l'],
    'cytochrome'     : ['1c2r_a','2mta_c'],
    'helical bundle' : ['1bbh_a','1bge_b','1rcb','1aep'],
    'ef-hand'        : ['1osa','2sas'],
    'other alpha'    : ['1hom','1lga_a','2hpd_a'],
    'tim barrel'     : ['1chr_a','2mnr','3rub_l'],
    'hydrolase'      : ['1crl','1tah_a'],
    'thieredoxin'    : ['1aba','1dsb_a','1gpl_a'],
    'ribonuclease'   : ['1atn_a','1hrh_a'],
    'open sheet'     : ['3chy','2ak3_a','1gky','2cmd','1eaf','2gbp','1mio_c','2pia','1gal','1npx'],
    'mixed'          : ['2hhm_a'],
    'small'          : ['1hip','1isu_a'],
    'ig'             : ['1fc1_a','2fbj_l'],
    'ig-like'        : ['1cid','1pfc','1ten','1tlk','3cd4','3hla_b'],
    'copredoxin'     : ['1aaj','2afn_a','2aza_a'],
    'virus'          : ['4sbv_a','1bbt_1'],
    'lectin-like'    : ['1sac_a'],
    'ob-fold'        : ['1lts_d'],
    'trefoil'        : ['1tie','8i1b'],
    'trypsin'        : ['1arb','2sga','2snv'],
    'lipocalin'      : ['1mdc','1mup'],
    'propeller'      : ['2sim'],
    'other beta'     : ['1cau_b','2omf'],
    'ub fold'        : ['1fxi_a'],
    'cystatin'       : ['1cew','1stf_i'],
    'sh2'            : ['2pna'],
   'other alpha+beta': ['2sar_a','1onc','5fd1']
}



# map id to name of class
FISCHER_ID_CLASS_DICT = {
    '1dxt_b' : 'alpha',
    '1cpc_l' : 'alpha',
    '1c2r_a' : 'alpha',
    '2mta_c' : 'alpha',
    '1bbh_a' : 'alpha',
    '1bge_b' : 'alpha',
    '1rcb'   : 'alpha',
    '1aep'   : 'alpha',
    '1osa'   : 'alpha',
    '2sas'   : 'alpha',
    '1hom'   : 'alpha',
    '1lga_a' : 'alpha',
    '2hpd_a' : 'alpha',
    '1chr_a' : 'alpha/beta',
    '2mnr'   : 'alpha/beta',
    '3rub_l' : 'alpha/beta',
    '1crl'   : 'alpha/beta',
    '1tah_a' : 'alpha/beta',
    '1aba'   : 'alpha/beta',
    '1dsb_a' : 'alpha/beta',
    '1gpl_a' : 'alpha/beta',
    '1atn_a' : 'alpha/beta',
    '1hrh_a' : 'alpha/beta',
    '3chy'   : 'alpha/beta',
    '2ak3_a' : 'alpha/beta',
    '1gky'   : 'alpha/beta',
    '2cmd'   : 'alpha/beta',
    '1eaf'   : 'alpha/beta',
    '2gbp'   : 'alpha/beta',
    '1mio_c' : 'alpha/beta',
    '2pia'   : 'alpha/beta',
    '1gal'   : 'alpha/beta',
    '1npx'   : 'alpha/beta',
    '2hhm_a' : 'other',
    '1hip'   : 'other',
    '1isu_a' : 'other',
    '1fc1_a' : 'beta',
    '2fbj_l' : 'beta',
    '1cid'   : 'beta',
    '1pfc'   : 'beta',
    '1ten'   : 'beta',
    '1tlk'   : 'beta',
    '3cd4'   : 'beta',
    '3hla_b' : 'beta',
    '1aaj'   : 'beta',
    '2afn_a' : 'beta',
    '2aza_a' : 'beta',
    '4sbv_a' : 'beta',
    '1bbt_1' : 'beta',
    '1sac_a' : 'beta',
    '1lts_d' : 'beta',
    '1tie'   : 'beta',
    '8i1b'   : 'beta',
    '1arb'   : 'beta',
    '2sga'   : 'beta',
    '2snv'   : 'beta',
    '1mdc'   : 'beta',
    '1mup'   : 'beta',
    '2sim'   : 'beta',
    '1cau_b' : 'beta',
    '2omf'   : 'beta',
    '1fxi_a' : 'alpha+beta',
    '1cew'   : 'alpha+beta',
    '1stf_i' : 'alpha+beta',
    '2pna'   : 'alpha+beta',
    '2sar_a' : 'alpha+beta',
    '1onc'   : 'alpha+beta',
    '5fd1'   : 'alpha+beta'
}

# map name of class to list of ids
FISCHER_CLASS_IDLIST_DICT = {
    'alpha'       : ['1dxt_b','1cpc_l','1c2r_a','2mta_c', '1bbh_a','1bge_b','1rcb','1aep','1osa','2sas', '1hom','1lga_a','2hpd_a'],
    'alpha/beta'  : ['1chr_a','2mnr','3rub_l','1crl','1tah_a','1aba','1dsb_a','1gpl_a', '1atn_a','1hrh_a','3chy','2ak3_a','1gky','2cmd','1eaf','2gbp','1mio_c','2pia','1gal','1npx'],
    'other'       : ['2hhm_a','1hip','1isu_a'],
    'beta'        : ['1fc1_a','2fbj_l', '1cid','1pfc','1ten','1tlk','3cd4','3hla_b', '1aaj','2afn_a','2aza_a','4sbv_a','1bbt_1', '1sac_a','1lts_d', '1tie','8i1b', '1arb','2sga','2snv', '1mdc','1mup', '2sim', '1cau_b','2omf'],
    'alpha+beta'  : ['1fxi_a', '1cew','1stf_i', '2pna', '2sar_a','1onc','5fd1']
}

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-rcv] "
                     " <query_id> <tabsearch.outputfile>\n")
    sys.stderr.write('  -c class level not fold level evaluation\n')
    sys.stderr.write('  -r higher scores are better (for MSVNS4MaxCMO)\n')
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.exit(1)

    
def main():
    """
    main for tsevalfischer.py

    Usage: tsevalfischer.py [-crv]  <query_id> <tabsearch.outputfile>

    
    -c evaluate at class level rather than default fold level
    -v turns on debug output to stderr
    -r higher scores are better (rather than default of more negative
       scores better). MSVNS4MaxCMO requires this (QP tableau search
       gives lower (more negative) scores for better matches)

    <query_id> is the PDB id (e.g. 1CRL or 1C2R_A) of the query structure
    
    <tabsearch.outputfile> is the output from the tabsearchqpml.file,
    which is a text file where each line is identifier then whitespace
    then score, sorted by score from most negative to least negative e.g.

    d1xksa_ -35.99999999
    d3sila_ -35.99999999
    ....
    d2mhua_ -0.499999999

    ie this means the best hit as the top of the file, and worst at bottom.
    May be specified as - for stdin.

    The table of positive and false negative rates is printed to stdout.
    """
    global verbose
    verbose = False
    use_class = False
    reverseflag = False
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "vrc?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        elif opt == "-r": # reversed: higher scores better
            reverseflag = True
        elif opt == '-c': # class not fold level evaluation
            use_class = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    query_id = args[0]
    tabsearch_file = args[1]

    if use_class:
        goldstd_ids = FISCHER_CLASS_IDLIST_DICT[FISCHER_ID_CLASS_DICT[query_id.lower()]]
    else:
        goldstd_ids = FISCHER_FOLD_IDLIST_DICT[FISCHER_ID_FOLD_DICT[query_id.lower()]]
    goldstd_ids = [pdbid.upper() for pdbid in goldstd_ids]

    if verbose:
        sys.stderr.write('parsing search results...\n')
    if tabsearch_file == '-':
        tabsearch_fh = sys.stdin
    else:
        tabsearch_fh = open(tabsearch_file)
    (searchresult,commentlist) = parse_searchresult(tabsearch_fh, reverseflag)
    if tabsearch_file != '-':
        tabsearch_fh.close()

    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('# results from:\n')
    for line in commentlist:
        sys.stdout.write('#  ')
        sys.stdout.write(line) # identifying information about search run
    sys.stdout.write('\n')
    eval_fn(goldstd_ids, searchresult, reverseflag)
    
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
