#!/usr/bin/env python

import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import QMCaller


caller = QMCaller (charge=-1, method="PM3", program="gaussian")
caller.Run ()
