#!/usr/bin/env python
###############################################################################
#
# star2auctable.py - convert StAR results.txt format to table of differences
#
# File:    star2auctable.py
# Author:  Alex Stivala
# Created: April 2010
#
# $Id: star2auctable.py 3681 2010-05-17 06:05:04Z alexs $
###############################################################################
"""
 star2auctable.py - convert StAR output to table of AUC differences

 Usage: star2auctable.py [-p pvalue] results.txt conf_intervals.txt reference-method

        -p pvalue - p-value to use for signficant differences, default 0.05

        results.txt is the output of StAR to read, a matrix of AUC differences
          (upper triangle) and p-values (lower triangle)

        conf_intervals.txt is the StAR conf_intervals.txt output file.

        reference-method is the method we want to measure the others
           relative to i.e we produce a table of AUC difference from
            this method with p-values
            (don't include the quotes that StAR always includes)

 Output is to stdout.
 
 The reference for StAR is:

 Vergara, Normabuena, Ferrada, Slater and Melo 'StAR: a simple tool for
 the statistical comparison of ROC curves' BMC Bioinformatics 2008 9:265

"""


import sys,os
import getopt

import numpy

def parse_results_file(results_fh):
    """
    Parse the StAR results.txt file, a matrix of AUC differences
    (upper triangle) and p-values (lower triangle).
    The file is tab-delimited

    Parameters:
       results_fh - open filehandle to read results.txt from

    Return value:
      tuple (resarray, methodlist)
        where resarray is the matrix parsed, with delta AUC in uppper
        and p-value in lower triange, n square for n methods.
        methodlist is list of method names, position in list is index
         in resarray for that method
    """
    for line in results_fh:
        sline = line.split('\t')
        if len(sline) < 2:
            continue
        if line[0] == '\t':
            methodlist = line.split('\t')[1:]
            n = len(methodlist)
            resarray = numpy.zeros((n,n))
            i = 0
            continue
        j = 0
        for v in sline[1:]:
            if i != j:
                resarray[i,j] = v
            j += 1
        i += 1

    # remove quotes from method names (and newline from last one)
    methodlist = [method.rstrip().lstrip('"').rstrip('"') for method in methodlist]
    return (resarray, methodlist)


def parse_conf_intervals(ci_fh):
    """
    Parse the StAR conf_intervals.txt file, each row a pair of methods
    with AUC difference (this time WITH sign, so we know which is better
    and which worse) and confidence interval

    Parameters:
        ci_fh - open filehandle to read conf_intervals.txt from

    Return value:
       dict { (method1,method2) : (auc_difference, cilower, ciupper) }

         mapping pair of methods to difference in AUC (method1 - method2),
           and lower and upper confidence interval value
    """
    ci_dict = {}
    lineno = 1
    for line in ci_fh:
        if lineno == 1:
            lineno += 1
            continue
        sline = line.split('\t')
        (method1,method2) = sline[0].split('/')
        method1 = method1.lstrip('"').rstrip('"')
        method2 = method2.lstrip('"').rstrip('"')
        deltaAUC = float(sline[1])
        cipair = sline[2] #  ( -0.0642863 , -0.0410837 )
        cilower = cipair.split(' ')[1]
        ciupper = cipair.split(' ')[3]
        ci_dict[(method1,method2)] = (deltaAUC, cilower, ciupper)
        lineno += 1
    return ci_dict
        
    

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage " + progname + "[-p pvalue] results.txt conf_intervals.txt reference-method\n")
    sys.exit(1)

def main():
    """
    main for star2auctable.py
    see usage message in header comment
    """
    sigpvalue = 0.05
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "p:h?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-p":
            sigpvalue = float(arg)
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 3:
        usage(os.path.basename(sys.argv[0]))

    resultsfilename = args[0]
    confintervalsfilename = args[1]
    referencemethod = args[2]

    (resarray, methodlist) = parse_results_file(open(resultsfilename))

    confintervals_dict = parse_conf_intervals(open(confintervalsfilename))


    q = 0
    for q in xrange(len(methodlist)):
        if methodlist[q] == referencemethod:
            j = q
            break
    if q >= len(methodlist)-1:
        sys.stderr.write("Method %s not found\n" % referencemethod)
        sys.exit (1 )

    notdiff_list = []
    for i in xrange(len(methodlist)):
        if methodlist[i] != referencemethod:
            try:
                ci_tuple = confintervals_dict[(referencemethod,methodlist[i])]
                signedDeltaAUC = ci_tuple[0]
            except KeyError:
                ci_tuple = confintervals_dict[(methodlist[i],referencemethod)]
                signedDeltaAUC = -ci_tuple[0]
            
            if i < j:
                deltaAUC = resarray[i,j]
                pvalue = resarray[j,i]
            else:
                pvalue = resarray[i,j]
                deltaAUC = resarray[j,i]

            assert(deltaAUC - abs(signedDeltaAUC) < 1e-5)
            
            if pvalue < sigpvalue:
                sys.stdout.write("%s\t%5.4f\t%5.4g\t%5.4f\n" % (methodlist[i], deltaAUC, pvalue, signedDeltaAUC))
            else:
                notdiff_list.append(methodlist[i])
    

    if len(notdiff_list) > 0:
#         sys.stdout.write("Not signficantly different from %s at p = %5.4f:\n" % (referencemethod, sigpvalue))
#         for m in notdiff_list:
#             sys.stdout.write("%s\n" % m)

        # output in same column format as different ones, to allow sort or
        # other easy parsing of output
        sys.stdout.write("%s\t%4.3f\t%5.4g\t%4.3f\n" %
                         (', '.join([referencemethod]+notdiff_list),
                          0, sigpvalue, 0))

        
if __name__ == "__main__":
    main()
