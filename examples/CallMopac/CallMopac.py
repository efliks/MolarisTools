#!/usr/bin/env python

import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCallerMopac


# . QMMM on
callerOn = QMCallerMopac (charge=-1, method="PM3", cosmo=False, qmmm=True, fileAtoms="mol.in", fileMopacError="run_with_qmmm.err", fileMopacInput="run_with_qmmm.mop", fileMopacOutput="run_with_qmmm.out", fileForces="forces_with_qmmm.out", fileTrajectory="qm_with_qmmm.xyz")
callerOn.Run ()

# . QMMM off
callerOff = QMCallerMopac (charge=-1, method="PM3", cosmo=True, qmmm=False, fileAtoms="mol.in")
callerOff.Run ()
