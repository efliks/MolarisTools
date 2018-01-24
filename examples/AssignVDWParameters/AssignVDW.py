#-------------------------------------------------------------------------------
# . File      : AssignVDW.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
"""Script reads a PDB file and assigns VDW parameters (A, B) to every atom."""

import sys, os

pathMolarisTools = os.path.join (os.environ["HOME"], "devel", "MolarisTools")
sys.path.append (pathMolarisTools)

from MolarisTools.Library  import AminoLibrary, ParametersLibrary
from MolarisTools.Parser  import PDBFile


pdb = PDBFile ("mbr_ts.pdb")
params = ParametersLibrary (os.path.join (pathMolarisTools, "data", "parm.lib"))
library = AminoLibrary (logging=True, verbose=False, unique=False, filename=os.path.join (pathMolarisTools, "data", "amino98_custom_small.lib"))

exclusions = ("HOH", "3NY", "SO4", )

for residue in pdb.residues:
    if not (residue.label in exclusions):
        resinlib = library[residue.label]
        for atom in residue.atoms:
            for libatom in resinlib.atoms:
                if (libatom.atomLabel == atom.label):
                    break
            for paramtype in params.vdws:
                if (paramtype.atomType == libatom.atomType):
                    break
            print ("%3s  %5d  %5s  %4d  %3s  %10.3f  %12.3f" % (residue.label, residue.serial, atom.label, atom.serial, libatom.atomType, paramtype.attractive, paramtype.repulsive))
