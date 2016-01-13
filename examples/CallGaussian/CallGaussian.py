#!/usr/bin/env python

import sys
sys.path.append ("/home/rcf-40/feliks/local/opt/MolarisTools")

from MolarisTools import QMCallerGaussian


caller = QMCallerGaussian (
    qmmm         = True           ,
    charge       = 0              ,
    method       = "B3LYP/6-31G*" ,
    fileAtoms    = "mol.in"       ,
    pathGaussian = "/usr/usc/gaussian/default/g03/g03" ,
        )
caller.Run ()
