#-------------------------------------------------------------------------------
# . File      : GenerateComponent.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os

from  MolarisTools.Scripts import AminoComponents_FromPDB
from  MolarisTools.Parser  import PDBFile


pdbfile = "avogadro.pdb"
components = AminoComponents_FromPDB (pdbfile)

component = components[0]
component.Write (showGroups=True, showLabels=True, sortGroups=True, filename="before.lib")


# . Calculate quantum charges for the component
pdb          = PDBFile (pdbfile)
geometry     = pdb.residues[0]
pathGaussian = os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03")

component.CalculateCharges (geometry, charge=0, ncpu=1, pathGaussian=pathGaussian)

component.Write (showGroups=True, showLabels=True, sortGroups=True, filename="after.lib")
