#!/usr/bin/env python

import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCallerGaussian


caller = QMCallerGaussian (charge=-1, method="PM3", fileAtoms="mol.in")
caller.Run ()
