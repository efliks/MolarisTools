#!/usr/bin/env python

import sys, os
sys.path.append ("/home/rcf-40/feliks/local/opt/MolarisTools")

from MolarisTools import QMCallerGaussian


# . Prepare the environment for Gaussian
env = os.environ
env["GAUSS_EXEDIR"] = "/usr/usc/gaussian/default/g03"
env["GAUSS_SCRDIR"] = "/tmp"

caller = QMCallerGaussian (
    env          = env            ,
    qmmm         = True           ,
    charge       = 0              ,
    method       = "B3LYP/6-31G*" ,
    fileAtoms    = "mol.in"       ,
    pathGaussian = "/usr/usc/gaussian/default/g03/g03" ,
        )
caller.Run ()
