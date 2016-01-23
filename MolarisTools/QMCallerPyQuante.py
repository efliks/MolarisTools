#-------------------------------------------------------------------------------
# . File      : QMCallerPyQuante.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import sys, os.path
sys.path.append (os.path.join (os.environ["HOME"], "local", "opt", "PyQuante-1.6.5", "build"))

from Units              import symbolToAtomicNumber
from QMCaller           import QMCaller

# . Generate file pyquante.log every time they are loaded (?)
from PyQuante           import SCF
from PyQuante.Molecule  import Molecule


class QMCallerPyQuante (QMCaller):
    """A class to provide communication between Molaris and PyQuante."""

    # . Options specific to PyQuante
    defaultAttributes = {
        "method" : "MINDO3" ,
        }
    defaultAttributes.update (QMCaller.defaultAttributes)

    def __init__ (self, **keywordArguments):
        super (QMCallerPyQuante, self).__init__ (**keywordArguments)
        # . Construct a PyQuante molecule
        atoms = []
        for atom in (self.molaris.qatoms + self.molaris.latoms):
            coordinates = (atom.x, atom.y, atom.z)
            pyqAtom     = (symbolToAtomicNumber[atom.symbol], coordinates)
            atoms.append (pyqAtom)
        self.molecule = Molecule ("molecule", atoms, units="Angstrom", charge=self.charge)


    def Run (self):
        calculation = SCF (self.molecule, method=self.method)
        calculation.iterate ()
        # . TODO


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
