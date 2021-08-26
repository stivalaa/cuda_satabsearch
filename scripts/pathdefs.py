###############################################################################
#
# pathdefs.py - Definitions of pathnames: edit for your system.
#
# File:    pathdefs.py
# Author:  Alex Stivala
# Created: August 2008
#
# $Id: pathdefs.py 3485 2010-03-17 04:48:13Z alexs $
# 
###############################################################################
"""
Locations of directory hierarchies for ASTRAL SCOP etc.
Edit appropriately for your system.

SCOP and ASTRAL data is obtained using the Bio.SCOP library (Casbon et
al 2006 'A high level interface to SCOP and ASTRAL implemented in
Python' BMC Bioinformatics 7:10) and depends on having the data
downloaded, in SCOP_DIR (defined below).

Downloaded SCOP files from

http://scop.mrc-lmb.cam.ac.uk/scop/parse/index.html

and ASTRAL files (in scopseq-1.73) from

http://astral.berkeley.edu/scopseq-1.73.html

The files downlaoded are:

/local/charikar/SCOP/:
dir.cla.scop.txt_1.73
dir.des.scop.txt_1.73
dir.hie.scop.txt_1.73

/local/charikar/SCOP/scopseq-1.73:
astral-scopdom-seqres-all-1.73.fa
astral-scopdom-seqres-sel-gs-bib-95-1.73.id

Other files there are indices built by Bio.SCOP when first used.

"""

# location of ASTRAL PDB-style coordinate files divided hierarchy
ASTRAL_ROOT = "/usr/local/ASTRAL/pdbstyle-1.75"

# location of SCOP dir files
#SCOP_DIR     = "/usr/local/SCOP"
SCOP_DIR     = "/home/alexs/SCOP"

# SCOP version to use
SCOP_VERSION = 1.75
