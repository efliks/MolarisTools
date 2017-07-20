#!/bin/bash

CUSTOM_PATH=./libs
export OUT_DIR=.

export EVB_LIB=${CUSTOM_PATH}/evb_poll_clean.lib
export AMINO_LIB=${CUSTOM_PATH}/amino98_custom_small.lib
export PARM_LIB=${CUSTOM_PATH}/parm.lib
export SOLVENT_OPT=${CUSTOM_PATH}/solvent.opt
export PDB_CONECT=${CUSTOM_PATH}/pdb_dictionary


for f in evb*inp ; do molaris_hpc9.15_apr06  $f > `basename $f .inp`.out  ; done
