#!/bin/bash
#===================================================
export OUT_DIR=.
CUSTOM_PATH=.

export EVB_LIB=${CUSTOM_PATH}/evb_poll_clean.lib
export AMINO_LIB=${CUSTOM_PATH}/amino98_custom_small.lib
export PARM_LIB=${CUSTOM_PATH}/parm.lib
export SOLVENT_OPT=${CUSTOM_PATH}/solvent.opt
export PDB_CONECT=${CUSTOM_PATH}/pdb_dictionary
#===================================================

for f in evb_heat_*inp ; do molaris_hpc9.15 $f >`basename $f .inp`.out ; done
