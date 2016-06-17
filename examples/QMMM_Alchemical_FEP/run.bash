#!/bin/bash

#===================================================
export OUT_DIR=.
CUSTOM_PATH=../../libraries

export EVB_LIB=${CUSTOM_PATH}/evb_poll_clean.lib
export AMINO_LIB=${CUSTOM_PATH}/amino98_custom_small.lib
export PARM_LIB=${CUSTOM_PATH}/parm.lib
export SOLVENT_OPT=${CUSTOM_PATH}/solvent.opt
export PDB_CONECT=${CUSTOM_PATH}/pdb_dictionary
#===================================================

f=fep.inp
molaris_hpc9.15 $f > `basename $f .inp`.out
