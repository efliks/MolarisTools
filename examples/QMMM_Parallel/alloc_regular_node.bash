#!/bin/bash
#SBATCH --job-name="qmmm"
#SBATCH --output="qmmm.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --export=ALL
#SBATCH -t 47:59:59


# . Modules and paths
module purge
module load  gaussian

export SCRDIR=/scratch/$USER/$SLURM_JOB_ID


#===================================================
export OUT_DIR=.
CUSTOM_PATH=${HOME}/DNA_polymerase/libs

export EVB_LIB=${CUSTOM_PATH}/evb_poll_clean.lib
export AMINO_LIB=${CUSTOM_PATH}/amino98_custom_small.lib
export PARM_LIB=${CUSTOM_PATH}/parm.lib
export SOLVENT_OPT=${CUSTOM_PATH}/solvent.opt
export PDB_CONECT=${CUSTOM_PATH}/pdb_dictionary
#===================================================


for f in evb_equil_*.inp ; do
    # . Copy files to scratch space
    cp $f                     $SCRDIR/
    cp call_gauss_lra_para.py $SCRDIR/
    cp sol_r.pdb              $SCRDIR/
    # . Prepare restart file
    resfile=`awk '/rest_in/ { print $2; }' $f`
    if [ "$resfile" ] ; then
        cp $resfile $SCRDIR/
        sed -i "/rest_in/ s/evb_equil_...\///" $SCRDIR/$f
    fi
    # . Run simulation
    cd $SCRDIR/
    molaris_hpc9.15_parallel_qmmm  $f  >  `basename $f .inp`.out
    cd -
    # . Copy results
    cat   $SCRDIR/qm.xyz   >> qm.xyz
    cp    $SCRDIR/job*.inp  .
    cp    $SCRDIR/job*.log  .
    cp    $SCRDIR/mol.in*   .
    cp    $SCRDIR/d.o*      .
    cp -r $SCRDIR/`basename $f .inp`      .
    cp    $SCRDIR/`basename $f .inp`.out  .
    rm -r $SCRDIR/*

    # . Check for a failed Gaussian job
    grep -q "Normal termination" job*.log ; if ! [ $? -eq 0 ] ; then break ; fi
    # . Check for a failed Molaris job
    grep -q "NORMAL TERMINATION" `basename $f .inp`.out ; if ! [ $? -eq 0 ] ; then break ; fi
done
