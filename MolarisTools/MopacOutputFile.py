#-------------------------------------------------------------------------------
# . File      : MopacOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine, WriteData
import collections

Atom     = collections.namedtuple ("Atom"     , "symbol x y z charge")
Force    = collections.namedtuple ("Force"    , "x y z")
ScanStep = collections.namedtuple ("ScanStep" , "Efinal atoms forces mcharges charges")


class MopacOutputFile (object):
    """A class to read a MOPAC output file."""

    def __init__ (self, filename="run.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _GetGradientLine (self, openfile):
        line     = next (openfile)
        tokens   = TokenizeLine (line, converters=[int, int, None, None, None, float, float, None])
        gradient = tokens[6]
        return gradient


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get the number of atoms
                if line.count ("TOTAL NO. OF ATOMS:"):
                    tokens = TokenizeLine (line, converters=[int], reverse=True)
                    natoms = tokens[0]

                # . Read gradients
                elif line.count ("FINAL  POINT  AND  DERIVATIVES"):
                    # . Skip the next two lines
                    next (lines)
                    next (lines)
                    self.forces = []
                    for i in range (natoms):
                        force = Force (x=self._GetGradientLine (lines) * GRADIENT_TO_FORCE, y=self._GetGradientLine (lines) * GRADIENT_TO_FORCE, z=self._GetGradientLine (lines) * GRADIENT_TO_FORCE)
                        self.forces.append (force)

                # . Get the final total energy
                # elif line.count ("TOTAL ENERGY"):
                #     tokens = TokenizeLine (line, converters=[None, None, None, float, None])
                #     self.Etotal = tokens[3] * EV_TO_KCAL_MOL

                # . Read the final QM energy
                elif line.count ("FINAL HEAT OF FORMATION"):
                    tokens = TokenizeLine (line, converters=[None, None, None, None, None, float, None, None, float, None])
                    self.Efinal = tokens [5]

                # . Read ESP (= Merz-Kollman) charges
                elif line.count ("ELECTROSTATIC POTENTIAL CHARGES"):
                    # . Skip the next two lines
                    next (lines)
                    next (lines)
                    self.mkcharges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.mkcharges.append (charge)

                # . Read Mulliken charges
                elif line.count ("MULLIKEN POPULATIONS AND CHARGES"):
                    # . Skip the next two lines
                    next (lines)
                    next (lines)
                    self.charges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float, float])
                        charge = tokens[3]
                        self.charges.append (charge)
        except StopIteration:
            pass


    def WriteMolarisForces (self, filename="forces.out", Eref=0., useESPCharges=False):
        """Write a file in the Molaris-suitable format."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
