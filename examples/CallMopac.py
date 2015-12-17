#!/usr/bin/env python

import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCallerMopac


caller = QMCallerMopac (charge=-1, method="PM3")
caller.Run ()
