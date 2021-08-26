#!/bin/sh
#
# File:    normalize_all.sh
# Author:  Alex Stivala
# Created: October 2008
#
# normalize_all.sh - build all normalized scores files
#
# Usage: normalize_all.sh 
#
# Run from the qptabsearch/ directory: processes files under
# results/fischer, results/nh3d, maxcmo_results/fishcer, maxcmo_results/nh3d
# etc., makeing norm1/ norm2/ norm3/ subdirectories under those,
# using also tableaux and contact map input files under data/
# via scripts normalize_tabmatch.py and normalize_msvns4maxcmo.py
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_all.sh 2092 2009-03-09 23:19:14Z astivala $

TABMATCH_FISCHER_DIR=results/fischer
TABMATCH_NH3D_DIR=results/nh3d
MSVNS_FISCHER_DIR=maxcmo_results/fischer
MSVNS_NH3D_DIR=maxcmo_results/nh3d

TABLEAUX_FISCHER_DIR=data/fischer_db
TABLEAUX_NH3D_DIR=data/nh3d
CM_FISCHER_DIR=data/fischer_cm
CM_NH3D_DIR=data/nh3d_cm

TABMATCH_QUERY200_DIR=results/query200
TABLEAUXDB=/local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.full.ascii

TABSEARCH_QUERY200_DIR=other_results/TableauSearch/query200
TABSEARCHDB=/local/charikar/TableauSearchDB

NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing QP tableau search Fischer results with norm ${norm}..."
    outdir=${TABMATCH_FISCHER_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_FISCHER_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
        normalize_tabmatch.py ${norm} ${qid} ${TABLEAUX_FISCHER_DIR} < ${infile} > ${outfile}
    done
    echo "normalizing QP tableau search Nh3D results with norm ${norm}..."
    outdir=${TABMATCH_NH3D_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_NH3D_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
        # dodgy: remove periods so that CATH id fits in 8 chars... 
        # luckily we get no duplicates...
        cathid=`echo ${qid} | tr -d .`
        normalize_tabmatch.py ${norm} ${cathid} ${TABLEAUX_NH3D_DIR} < ${infile} > ${outfile}
    done

    echo "normalizing MSVNS Fischer results with norm ${norm}..."
    outdir=${MSVNS_FISCHER_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${MSVNS_FISCHER_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
        normalize_msvns4maxcmo.py ${norm} ${qid} ${CM_FISCHER_DIR} < ${infile} > ${outfile}
    done
    echo "normalizing MSVNS Nh3D results with norm ${norm}..."
    outdir=${MSVNS_NH3D_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${MSVNS_NH3D_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
        # dodgy: remove periods so that CATH id fits in 8 chars... 
        # luckily we get no duplicates...
        cathid=`echo ${qid} | tr -d .`
        normalize_msvns4maxcmo.py ${norm} ${cathid} ${CM_NH3D_DIR} < ${infile} > ${outfile}
    done
    echo "normalizing QP tableau search query200 results with norm ${norm}..."
    outdir=${TABMATCH_QUERY200_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_QUERY200_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
#        echo $qid
        normalize_tabmatch.py ${norm} ${qid} ${TABLEAUXDB} < ${infile} > ${outfile}
    done
    echo "normalizing TableauSearch query200 results with norm ${norm}..."
    outdir=${TABSEARCH_QUERY200_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABSEARCH_QUERY200_DIR}/*.scores
    do
        qid=`basename ${infile} .scores`
        outfile=${outdir}/${qid}.tabsearch.out
        tableausearchout2col.py < ${infile} | normalize_tabmatch.py ${norm} ${qid} ${TABSEARCHDB} > ${outfile}
    done
done

