#!/usr/bin/env python
#-------------------------------------------------------------------------------
# . File      : caller_gaussian.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import sys, os

# . Define a path to the MolarisTools library or comment this out and 
# . place MolarisTools in the current directory
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))


from MolarisTools import QMCallerGaussian

# . Modify the number of CPUs, memory requirements (usually 1 GB 
# . per CPU is fine), electronic state and method
#
# . qmmm=True  for electrostatic embedding (requires use_qm_reg2w_force in the Molaris input file)
#
# . qmmm=False for mechanical    embedding (requires cl_elec in the Molaris input file)
callerGaussian = QMCallerGaussian (
        ncpu            =  1      ,
        memory          =  1      ,
        archive         =  True   ,
        restart         =  True   ,
        qmmm            =  True   ,
        charge          =    0             ,
        multiplicity    =    1             ,
        method          =  "B3LYP/6-31G*"  ,
        pathGaussian           = os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03") ,
        fileGaussianCheckpoint = os.path.join (os.environ["PWD"], "job.chk") ,
            )
callerGaussian.Run ()
