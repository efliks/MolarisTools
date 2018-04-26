#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : CallGaussian.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os

from MolarisTools.QMMM  import QMCallerGaussian


# . Modify the number of CPUs, memory requirements (usually 1 GB 
# . per CPU is fine), electronic state and method.
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
