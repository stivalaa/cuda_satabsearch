#!/usr/bin/env python
###############################################################################
#
# mkroctabs.py - Make table of AUC scores or slrtab from multiquery output file
#
# File:    mkroctabs.py
# Author:  Alex Stivala
# Created: February 2010
#
# $Id: mkroctabs.py 3625 2010-05-07 07:54:39Z alexs $
# 
###############################################################################

"""
Build a table of AUC scores or .slrtab for R from a multiquery
output file from a db
search program (in the two-column (dbid score) format), with a comment
('#' in first column) in this format:

# QUERYID = querid

or

# QUERY ID = querid

e.g.

# QUERY ID = d1ubia_

to delimit results for each separate query. The blastout2col.sh etc.
scripts create this format, as does cudaSaTabsearch.

This way the SCOP database only has to be loaded once, rather than
running tsevalfn.py or similar once for each query - loading the
database is the slowest part so this is much more efficient.

See usage in docstring for main()

"""

import sys,os
import getopt
from itertools import groupby
 
from Bio.SCOP import *

from pathdefs import SCOP_DIR,SCOP_VERSION
from tsevalutils import (iter_searchresult, get_goldstd_domains,
                         tsevalutils_set_verbose,is_true_positive,
                         get_betagrasp_containing_domains,
                         filter_domains_astral_nrpercent)


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def compute_auc_mann_whitney(searchresult, goldstd_domains,
               use_nonredundant, 
               scop, astral):
    """
    Compute the AUC by by means of the Mann-Whitney U statistic
    (Wilcoxon rank-sum test) divided by product of number of false
    positives and the total number of true positives possible
    in the data.

    The AUC score is the area under the ROC curve, This is computed by
    processing scores (from best to worst) marking them as true
    postivies (same fold/superfamily as query) or false positives
    (different fold/superfamily from query).  For each of the
    positives, the number of true positives with a better simliarity
    score is calculatd, and this is divided by the number of false
    postives. This is then divided by the total number of possible
    true postivies in the database (i.e. the number of mumbers in the
    SCOP fold/superfamily).

    
    Parameters:
       searchresult -list of (score, domainid) tuples parsed from pogram to eval
                     This is assumed to be sorted from best to worst score
       goldstd_domains - list of gold standard domains to eval against
       use_nonredundant - Bool if True use ASTRAL 95% nr subset
       scop - previously built Bio.SCOP Scop instance
       astral - previously build Bio.SCOP Astral instance

    Return value:
       The AUC for the searchresults evaluated against the 
       goldstd_domains.
    """
    # convert goldstd list of domainids to dictionary keyed by domainid
    # for fast lookup as we iterate through search results.
    # The dictionary is { domainid : True } (we don't have a value,
    # just need to quickly test for presence of domainid in gold std
    # postive list).
    goldstd_pos_dict = dict([(sid, True) for
                       sid in [scopdom.sid for scopdom in goldstd_domains]])

    total_tp_possible = len(goldstd_domains) 
#    print 'bbb total tp possible',total_tp_possible

    labelled_ranks = [] # list of (sid, score, rank, is_pos)
    i = 0
    fp_count = 0
    while i < len(searchresult):
        try:
            (score, domainid) = searchresult[i]
        except IndexError:
            sys.stderr.write('WARNNG: only %d search results, only got %d false positives\n' % (len(searchresult), fp_count))
            break
#        print 'zzz',score,domainid

        istruepos = is_true_positive(domainid, goldstd_pos_dict)
        if not istruepos:
            fp_count += 1
#            if fp_count == 1:
#                print 'eee first fp',domainid,'at ',i
        labelled_ranks.append((domainid, score, i, istruepos))
        i += 1

    if fp_count == 0:
        # if no false positives found, then set AUC to 1.0 as it means
        # the results are 'perfect' up to the number of scores reported.
        sys.stderr.write('WARNING: found zero false positives, set AUC=1.0\n')
        return 1.0

    tp_count = len([sid for (sid,score,rank,istruepos) in labelled_ranks 
                    if istruepos])

#    print 'ddd fp_count',fp_count
#    print 'fff tp_count',tp_count
    assert(tp_count == len(labelled_ranks) - fp_count)
    
#    for tup in labelled_ranks:
#        print 'xxx',tup


#     ########################################################################
#     # naive 'direct' computation 
#     total = 0
#     for i in xrange(len(labelled_ranks)):
#         count = 0
#         if not labelled_ranks[i][3]:
#             for j in xrange(i):
#                 if labelled_ranks[j][3]:
#                     count += 1
#             total += count
#     print 'nnn total = ',total
#     print 'ooo ', float(total) / ( float(fp_count) * float(total_tp_possible) )
#     #########################################################################

    fp_rank_sum = sum([rank+1 for (sid,score,rank,istruepos)
                       in labelled_ranks if not istruepos])
    mann_whitney_U = fp_rank_sum - fp_count*(fp_count+1)/2
#    print 'uuuu U =',mann_whitney_U

    tp_rank_sum = len(labelled_ranks)*(len(labelled_ranks)+1)/2 - fp_rank_sum
    mann_whitney_U2 = tp_rank_sum - tp_count*(tp_count+1)/2
#    print 'vvvv U2 = ',mann_whitney_U2

    assert(mann_whitney_U + mann_whitney_U2 == fp_count * tp_count)

#    auc = float(mann_whitney_U) / float(fp_count*total_tp_possible)
##    sys.stderr.write('xxx ' + str(fp_count) + ' ' + str(tp_count)+'\n')
    if (tp_count == 0):
        if (total_tp_possible == 0):
            sys.stderr.write('WARNING: found zero true positives, but 0 true positives exist anyway, set AUC=1.0\n')
            auc = 1.0
        else:
            sys.stderr.write('WARNING: found zero true positives, (but %s exist), set AUC=0.0\n' % total_tp_possible)
            auc = 0.0
    else:
        auc = float(mann_whitney_U) / float(fp_count *  tp_count)


#     #########################################################################
#     # average precision calculation
#     total = 0
#     for i in xrange(len(labelled_ranks)):
#         higher_tp_count = 1
#         if labelled_ranks[i][3]:
#             for j in xrange(i):
#                 if labelled_ranks[j][3]:
#                     higher_tp_count += 1
#             pos_count = i+1
#             prec = float(higher_tp_count) / float(pos_count)
#             print 'qqqq i prec = ',i,prec
#             total += prec
#     avgprec = total / float(len(labelled_ranks)) # paper says this; seems wrong
# #    avgprec = total / float(tp_count)
#     print 'pppp avgprec = ',avgprec
#     #########################################################################


    return auc


def write_slrtab(searchresult, goldstd_domains,
                 use_nonredundant, 
                 scop, astral):
    """
    Write the slrtab output to stdout. This is two columms with score in
    first column and label (0 or 1) in second column, for use with R
    ROCR package.
    
    Parameters:
       searchresult -list of (score, domainid) tuples parsed from pogram to eval
                     This is assumed to be sorted from best to worst score
       goldstd_domains - list of gold standard domains to eval against
       use_nonredundant - Bool if True use ASTRAL 95% nr subset
       scop - previously built Bio.SCOP Scop instance
       astral - previously build Bio.SCOP Astral instance

    Return value:
       None.
    """
    # convert goldstd list of domainids to dictionary keyed by domainid
    # for fast lookup as we iterate through search results.
    # The dictionary is { domainid : True } (we don't have a value,
    # just need to quickly test for presence of domainid in gold std
    # postive list).
    goldstd_pos_dict = dict([(sid, True) for
                       sid in [scopdom.sid for scopdom in goldstd_domains]])

    for (score, domainid) in searchresult:
        if goldstd_pos_dict.has_key(domainid):
            label = 1
        else:
            label = 0
        sys.stdout.write('%20.8f  %d\n' % (score,label))



#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname +" [-v] [-e version] [-p percent] [-n] [-l] [-b] [a|u] [-z score] [-s]\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -u evaluate at superfamily rather than fold level\n')
    sys.stderr.write('  -a evaluate at family rather than fold level\n')
    sys.stderr.write('  -e version: SCOP version (default %s)\n' %  SCOP_VERSION)
    sys.stderr.write('  -p percent: ASTRAL sequence percent nonredundant subset to use\n')
    sys.stderr.write('  -n negate scores\n')
    sys.stderr.write('  -l write table of values and labels intead of AUC\n')
    sys.stderr.write('  -b use the ProSMoS data set for validation of beta-grasp query\n')
    sys.stderr.write('  -z score : assign identifiers not present in the output a score of score\n')
    sys.stderr.write('  -s : skip self-query\n')
    sys.exit(1)

    
def main():
    """
    main for mkroctabs.py

    Usage: mkroctabs.py [-v] [-e version] [-p percent] [-n] [-l] [a|u] [-b] [-z score] [-s]

    
    -v turns on debug output to stderr

    -u use domains from the same superfamily as the query rather than from
       the same fold as the query as the gold standard.

    -a use domains from same family as the query rather tahn from the
       same fold as the query as the gold standard.

    -e version: SCOP version to use (default 1.73)

    -p percent: ASTRAL sequence nonredundant subset percentage to use
                (e.g. 95). If not specified, use all domains.

    -n negate all scores (so that most -ve score becomes highest score)

    -l instead of computing AUC, just
       write table with one column of values (scores) and other of 
       corresponding true labels (0/1 for same/different fold as query)
       for use with the R ROCR packge.

    -b use the PRoSMoS data set (Shi et al 2007) data set for validation
       of the beta-grasp substructure query

    -z score : any identifiers that are not present in the output but
               are in the gold standard data are given the specified score.
               This is for programs that do not assign a score to all
               domains, but only those above some threshold or the top n,
               or just cannot assign a score to some domains for some reason.
               This would generally be specified as some score lower 
               than all other scores.

    -s omit score of query against itself

    input on stdin is the two column (identifier score) output of the
    program to evaluate for that query as described in module header
    docstring.

    output on stdout is the table of AUC values (two column (queryid AUC)) 
    for all queries concatenated together, with header for use in R
    read.table(header=TRUE).
    """
    global verbose
    verbose = False
    
    negateflag = False
    use_prosmos_dataset = False
    use_nonredundant = False
    use_superfamily = False
    use_family = False
    bottom_score = None
    logflag = False
    nrpercent = None
    scop_version = SCOP_VERSION
    do_slrtab = False
    use_prosmos_dataset = False
    bottom_score = None
    skip_self_query = False

    try:
        opts,args = getopt.getopt(sys.argv[1:], "bauve:p:nlz:s?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
            tsevalutils_set_verbose(True) # tsevalutils.py verbose flag
        elif opt == "-u": # evaluate at superfamily not fold level
            use_superfamily = True
        elif opt == "-a": # evaluate at family not fold level
            use_family = True
        elif opt == "-p": # sequence identity percentage for nr subset
            nrpercent = int(arg)
            use_nonredundant = True
        elif opt == "-e": # SCOP version number
            scop_version = arg
        elif opt == "-n": # negate scores
            negateflag = True
        elif opt == "-l": # write slrtab instead of computing AUC
            do_slrtab = True
        elif opt == "-b": # use ProSMoS beta-grasp dataset
            use_prosmos_dataset = True
        elif opt == "-z": # score to give to domains that have no score
            bottom_score = float(arg)
        elif opt == "-s": # omit result of query struct against itself
            skip_self_query = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))

    if use_family and use_superfamily:
        sys.stderr.write('-a (family) and -u (superfamily) are mutually exclusive\n')
        usage(sys.argv[0])

    level = 'fold'
    if use_family:
        level = 'family'
    elif use_superfamily:
        level = 'superfamily'

    # read SCOP and ASTRAL data
    if verbose:
        sys.stderr.write('Reading SCOP data (version %s)...\n' % scop_version)
    scop = Scop(dir_path=SCOP_DIR,version=scop_version)
    astral = Astral(dir_path=SCOP_DIR,version=scop_version,scop=scop)


    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    if do_slrtab:
        sys.stdout.write('score  label\n')
        sys.stdout.write('#-----------\n')
    else:
        sys.stdout.write('queryid          AUC\n')
        sys.stdout.write('#-------------------\n')

    

    # get list of iterables each for same queryid.
    # iter_searchresult() is isterable of tuples (queryid, score, domainid)
    # groupby() requires that the iterable already has identical consecutive
    # queryids (first element of tuple) - iter_searchresult() should yield this
    # and does, when there is only one instance of the QUERY ID for each query,
    # but new version of cudaSaTabsearch does all queries in small structure
    # db, then in large structure db, so two QUERY ID for each query in the
    # output, so we sort by queryid first before groupby.
    query_group_iter = groupby(sorted(
        iter_searchresult(sys.stdin, multiquery=True,
                          skip_self_query=skip_self_query,
                          negateflag=negateflag) ),
                               lambda t : t[0])
    
    total = 0
    num = 0
    error_count = 0
    for (queryid, result_iter) in query_group_iter:
        if verbose:
            sys.stderr.write('processing query id %s\n' % queryid)
        if use_prosmos_dataset:
            goldstd_domains = get_betagrasp_containing_domains(scop)
            if verbose:
                sys.stderr.write("found %d betagrasp containing domains\n"
                                 % len(goldstd_domains))
            if use_nonredundant:
                goldstd_domains = filter_domains_astral_nrpercent(
                    goldstd_domains, scop, astral, nrpercent)
                if verbose:
                    sys.stderr.write("reduced to %d betagrasp containing domains in %d %% nonredundant subset\n" % (len(goldstd_domains), nrpercent))
        else:
            try:
                goldstd_domains = get_goldstd_domains(queryid.lower(), 
                                                      level,
                                                      use_nonredundant,
                                                      scop, astral,
                                                      nrpercent)
            except AttributeError:
                # AttributeError: 'NoneType' object has no attribute 'getAscendent'
                # when we get a None from looking up query domain identifiAer
                sys.stderr.write('ERROR: query sid %s not found in SCOP, skipping\n' % queryid.lower())
                error_count += 1
                continue
        
        allscores = [(score,dbid) for (queryid,score,dbid) in result_iter]

        if bottom_score != None:
            lowscore_domains = 0
            # build dictinoary for fast lookup of sids in search result
            searchresult_dict = dict([(sid.lower(), True) for (score,sid) in allscores])
            all_domains = scop.getRoot().getDescendents('domain')
            if use_nonredundant:
                all_domains = filter_domains_astral_nrpercent(all_domains,
                                                              scop, astral,
                                                              nrpercent)
            for scopdom in all_domains:
                if skip_self_query and scopdom.sid.lower() == queryid.lower():
                    continue # skip self query
                if not searchresult_dict.has_key(scopdom.sid.lower()):
                    lowscore_domains += 1
                    allscores.append((bottom_score,scopdom.sid.lower()))
            if verbose and lowscore_domains > 0:
                sys.stderr.write("(queryid %s): set score to %f for %d domains\n" % (queryid, bottom_score, lowscore_domains))

        if skip_self_query:
            # skipping query against itelf, need to remove it from gold standrd
            goldstd_domains = [d for d in goldstd_domains
                               if d.sid.lower() != queryid.lower()]
            
        if do_slrtab:
            sys.stdout.write('# QUERY ID = ' + queryid + '\n')
            auc = write_slrtab(sorted(allscores,reverse=True),
                               goldstd_domains,
                               use_nonredundant, scop, astral)
        else:
            auc = compute_auc_mann_whitney(
                sorted(allscores,reverse=True),
                goldstd_domains,
                use_nonredundant, scop, astral)
            total += auc
            num += 1
            sys.stdout.write('%s    %4.3f\n'  % (queryid.lower(), auc))

            avg = total / num

    if not do_slrtab:
        sys.stdout.write('# AVERAGE = %f\n' % (avg))
            
            
if __name__ == "__main__":
    main()
