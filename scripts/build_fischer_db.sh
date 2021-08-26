#!/bin/sh
#
# File:    build_fischer_db.sh
# Author:  Alex Stivala
# Created: September 2008
#
# build_fischer_db.sh - build tableaux database for Fischer data set
#
# Usage: build_fischer_db.sh outdir 
#
#   outdir is name of diretory which is created, and each tableau 
#   in ASCII format for use with tsrchd_sparse etc. is
#   created as a separate file in that directory, in format for input
#   for use with qptabmatch_allpairs.py for example
#
#   To stdout is written the ASCII format db of all the tableaux+dist matrices
#   (just all the ones written to outdir concatenated together with
#   blank line between each). 
#  
# builds database of tableaux, using pytableaycreate.py,
# for the Fischer data set (Fischer et al 1996 Pac. Symp. Biocomput. 300-318))
# This allows all-against-all (including redundant, so for n (=68)
# there are n*n (=4624) total comparions) with e.g. qptabmatch_allall.py,
#
#

# root of divided PDB hierarchy
PDBROOT=/local/charikar/pdb/pdb


# List of probe PDB ids from Fischer 1996 Table I
# Note several PDB ids obsoleted, so change to the replacments
FISCHER_S="1mdc 1mup 1npx 1cpc_l 1onc 2ak3_a 1osa 1atn_a 1pfc 1arb 2cmd 2pia 2pna 3rub_l 1bbh_a 2sar_a 1c2r_a 3cd4 1chr_a 1aep 1dxt_b 2mnr 2fbj_l 1lts_d 1gky 2gbp 1hip 1bbt_1 2sas 2mta_c 1fc1_a 1tah_a 2hpd_a 1rcb 1aba 1sac_a 1eaf 1dsb_a 2sga 1stf_i 2hhm_a 2afn_a 1aaj 1fxi_a 5fd1 1bge_b 1isu_a 3hla_b 1gal 3chy 1cau_b 2aza_a 1hom 1cew 1tlk 1cid 2omf 1crl 1lga_a 2sim 1mio_c 1ten 4sbv_a 1tie 8i1b 2snv 1hrh_a 1gp1_a"


# List of target fold PDB ids from Fischer 1996 Table I
# Note several PDB ids obsoleted, so change to the replacments
# this list corresponds to FISCHER_S ie FISCHER_P[i] is the target fold
# for probe FISCHER_S[i] for 0 < i < 67
FISCHER_P="1ifc 1rbp 3grs 1col_a 7rsa 1gky 4cpv 1atr 3hla_b 5ptp 6ldh 1fnb 1sha_a 6xia 2ccy_a 9rnt 1ycc 2rhe 2mnr 256b_a 1hbg 4enl 8fab_b 1bov_a 3adk 2liv 2hip_a 2plv1 2scp_a 1ycc 2fb4_h 1tca 2cpp 2gmf_a 1ego 2ayh 4cla 2trx_a 5ptp 1mol_a 1fbp_a 1aoz_a 1paz 1ubq 1iqz 2gmf_a 2hip_a 2rhe 3cox 2fox 1cau_a 1paz 1lfb 1mol_a 2rhe 2rhe 2por 1ede 2cyp 1nsb_a 2min_b 3hhr_b 2tbv_a 4fgf 4fgf 5ptp 1rnh 2trx_a"


# List of 68 probe sequences from Fischer 1996 Table II
# Note several PDB ids obsoleted, so change to the replacments
FISCHER_LIST="1dxt_b 1cpc_l 1c2r_a 2mta_c 1bbh_a 1bge_b 1rcb 1aep 1osa 2sas 1hom 1lga_a 2hpd_a 1chr_a 2mnr 3rub_l 1crl 1tah_a 1aba 1dsb_a 1gpl_a 1atn_a 1hrh_a 3chy 2ak3_a 1gky 2cmd 1eaf 2gbp 1mio_c 2pia 1gal 1npx 2hhm_a 1hip 1isu_a 1fc1_a 2fbj_l 1cid 1pfc 1ten 1tlk 3cd4 3hla_b 1aaj 2afn_a 2aza_a 4sbv_a 1bbt_1 1sac_a 1lts_d 1tie 8i1b 1arb 2sga 2snv 1mdc 1mup 2sim 1cau_b 2omf 1fxi_a 1cew 1stf_i 2pna 2sar_a 1onc 5fd1"

if [ $# -ne 1 ]; then
    echo "Usage: $0  outdir" 2>&1
    exit 1
fi
outdir=$1

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

# pytableaucreate.py options
tabopts="-35 -f -t dssp -p none"

first=1
for i in $FISCHER_LIST
do
    pdb=`echo $i | tr A-Z a-z`
    if [ `expr index $pdb _` -ne 0 ]; then
        # get chainid from e.g. 1BYO_B
        chainid=`expr substr $pdb 6 1`
        chainopt="-c $chainid"
        pdbid=`expr substr $pdb 1 4`_${chainid}
    else
        chainopt=""
        pdbid=`expr substr $pdb 1 4`
    fi
    pdb=`expr substr $pdb 1 4`
    div=`expr substr $pdb 2 2`
    pdbfile=${PDBROOT}/${div}/pdb${pdb}.ent.gz
    if [ $first -eq 0 ]; then
        echo
    else
        first=0
    fi
    pytableaucreate.py ${tabopts} ${chainopt} ${pdbfile} | tee ${outdir}/${pdbid}.tableaudistmatrix 
    # append distance matrix, removing identifier on first line
    pytableaucreate.py -d ${tabopts} ${chainopt} ${pdbfile} | awk 'NR > 1'| tee -a ${outdir}/${pdbid}.tableaudistmatrix 
done

