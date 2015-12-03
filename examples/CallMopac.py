#!/usr/bin/env python

import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import MolarisAtomsFile, MopacOutputFile


# . Read atoms.inp from Molaris
molaris = MolarisAtomsFile (filename="atoms.inp")

# . Write geometries for the purpose of viewing
molaris.WriteQM (link=True, append=True)

# . Run a MOPAC calculation
molaris.WriteMopacInput (filename="run.mop", method="PM3", charge=-1, cosmo=True)
mopacError  = open ("run.err", "w")
subprocess.check_call ([os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe"), "run.mop"], stdout=mopacError, stderr=mopacError)
mopacError.close ()

# . Convert the output file from MOPAC to forces.out
mopac = MopacOutputFile (filename="run.out")
mopac.WriteMolarisForces (filename="forces.out")
