#!/usr/bin/env python

import sys, os
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools import PDBFile, AminoComponents_FromPDB


# . Generate component from a PDB file
components = AminoComponents_FromPDB ("mlt.pdb")
component  = components[0]

component.Write (showGroups=True, showLabels=True, sortGroups=True, filename="before.lib")


# . Calculate quantum charges for the component
pdb          = PDBFile ("mlt.pdb")
geometry     = pdb.residues[0]
pathGaussian = os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03")

component.CalculateCharges (geometry, charge=-2, ncpu=2, pathGaussian=pathGaussian)

component.Write (showGroups=True, showLabels=True, sortGroups=True, filename="after.lib")
