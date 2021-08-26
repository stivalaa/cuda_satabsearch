#!/bin/bash
#
# File:    qptabmatch_pbs_script.sh
# Author:  Alex Stivala
# Created: July 2009
#
# PBS script for submitting QP tableau search jobs on tango.vpac.org
# requires PATH and PYTHONPATH already set up in environment 
#
# $Id: qptabmatch_pbs_script.sh 2716 2009-07-29 02:07:14Z astivala $
 	
#PBS -N QP_tableau_search

#PBS -l walltime=1:0:0

#PBS -v MKL_NUM_THREADS=1

module load python

cd $PBS_O_WORKDIR
set CONV_RSH = ssh


time qptabmatchstructs.sh /home/alexs/share/ASTRAL/pdbstyle-sel-gs-bib-95-1.75/ql/d1qlpa_.ent /home/alexs/share/ASTRAL/pdbstyle-sel-gs-bib-95-1.75/ql/d1qlpa_.ent


