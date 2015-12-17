#!/usr/bin/env python
#
# . A simple script to compare PM3 jobs in Mopac and Gaussian
#
import sys, os, subprocess
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import MopacOutputFile, GaussianOutputFile


# . Create a Gaussian input file from a Mopac input file
lines  = open ("run.mop").readlines ()
mol    = []
for line in lines[3:]:
    tokens = line.split ()
    symbol, x, y, z = tokens[0], float (tokens[1]), float (tokens[3]), float (tokens[5])
    entry  = (symbol, x, y, z)
    mol.append (entry)

natoms   = len (mol)
gaussian = ["# PM3 FORCE\n", "\n", "...\n", "\n", "-1 1\n"]
for (symbol, x, y, z) in mol:
    gaussian.append ("%2s    %8.4f    %8.4f    %8.4f\n" % (symbol, x, y, z))
gaussian.append ("\n")

output = open ("gauss_test.inp", "w")
output.writelines (gaussian)
output.close ()

# . Run Gaussian
pathGaussian = os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03")
subprocess.check_call ([pathGaussian, "gauss_test.inp"])

# . Read output files
mopac    = MopacOutputFile    ("run.out")
gaussian = GaussianOutputFile ("gauss_test.log")

# . Compare energies
print ("Efinal(mopac) = %8.4f       Efinal(gaussian) = %8.4f" % (mopac.Efinal, gaussian.Efinal))
print ("")

# . Compare forces
# . Note: MopacOutputFile is missing "atoms", so implement it.
for atom, forceMopac, forceGaussian in zip (gaussian.atoms, mopac.forces, gaussian.forces):
    print ("%2s    Force(mopac) = %8.4f %8.4f %8.4f     Force(gaussian) = %8.4f %8.4f %8.4f     Force(difference) = %8.4f %8.4f %8.4f" % \
    (atom.symbol, \
    forceMopac.x, forceMopac.y, forceMopac.z, \
    forceGaussian.x, forceGaussian.y, forceGaussian.z, \
    forceMopac.x - forceGaussian.x, forceMopac.y - forceGaussian.y, forceMopac.z - forceGaussian.z ))
print ("")

# . Compare charges
for atom, chargeMopac, chargeGaussian in zip (gaussian.atoms, mopac.charges, gaussian.charges):
    print ("%2s    Charge(mopac) = %8.4f     Charge(gaussian) = %8.4f     Charge(difference) = %8.4f" % (atom.symbol, chargeMopac, chargeGaussian, chargeMopac - chargeGaussian))
