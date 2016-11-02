#-------------------------------------------------------------------------------
# . File      : QChemOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine, WriteData
import collections, exceptions


Atom     = collections.namedtuple ("Atom"     , "symbol  x  y  z")
Force    = collections.namedtuple ("Force"    , "x  y  z")

_ATOMS_PER_LINE = 6


class EfieldFile (object):
    """A class to read an efield file from Q-Chem."""

    def __init__ (self, filename="efield.dat"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        field = []
        try:
            while True:
                line         = next (lines)
                (ex, ey, ez) = TokenizeLine (line, converters=[float, float, float])
                # . There is not need to multiply the gradients by -1, since Molaris does it after reading the d.o file.
                vector       = (
                    ex  *  HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  , 
                    ey  *  HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  , 
                    ez  *  HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  , )
                field.append (vector)
        except StopIteration:
            pass
        # . Close file
        lines.close ()
        self.field = field


#-------------------------------------------------------------------------------
class QChemOutputFile (object):
    """A class to read a Q-Chem output file."""

    def __init__ (self, filename="job.log"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get geometry
                #             Standard Nuclear Orientation (Angstroms)
                #    I     Atom           X                Y                Z
                # ----------------------------------------------------------------
                #    1      C       0.3610741620    -1.0892812370     0.3333499280
                #    2      H      -0.5915297190    -0.4548440630     0.3617850530
                # (...)
                # ----------------------------------------------------------------
                if   line.count ("Standard Nuclear Orientation (Angstroms)"):
                    for i in range (2):
                        next (lines)
                    atoms = []
                    line  = next (lines)
                    while not line.count ("----"):
                        (serial, symbol, x, y, z) = TokenizeLine (line, converters=[int, None, float, float, float])
                        atom   = Atom (symbol=symbol, x=x, y=y, z=z)
                        atoms.append (atom)
                        line   = next (lines)
                    self.atoms = atoms


                # . Get Mulliken charges
                #          Ground-State Mulliken Net Atomic Charges
                #
                #     Atom                 Charge (a.u.)
                #  ----------------------------------------
                #      1 C                    -0.564153
                #      2 H                     0.296853
                # (...)
                #  ----------------------------------------
                elif line.count ("Ground-State Mulliken Net Atomic Charges"):
                    for i in range (3):
                        next (lines)
                    charges = []
                    line    = next (lines)
                    while not line.count ("----"):
                        tokens = TokenizeLine (line, converters=[int, None, float])
                        charges.append (tokens[2])
                        line   = next (lines)
                    self.charges = charges


                # . Get self energy of point charges
                #
                # Charge-charge energy     =  -250.9297020579 hartrees
                elif line.startswith (" Charge-charge energy"):
                    tokens      = TokenizeLine (line, converters=[None, float], reverse=True)
                    self.Echrg  = tokens[1] * HARTREE_TO_KCAL_MOL


                # . Get final SCF energy
                #
                # SCF   energy in the final basis set = -1211.5551873917
                elif line.startswith (" SCF   energy in the final basis set"):
                    tokens      = TokenizeLine (line, converters=[float, ], reverse=True)
                    self.Efinal = tokens[0] * HARTREE_TO_KCAL_MOL


                # . Get gradients on QM atoms
                #
                # Calculating analytic gradient of the SCF energy
                # Gradient of AOints
                #            1           2           3           4           5           6
                #    1   5.1830948   0.7198311  -1.3316471   5.4310464   8.1701661 -15.4180387
                #    2  -1.7031303  -3.3832442   7.5553033  -3.1040094  -0.7271863  -5.9138090
                #    3  -6.0741297  -1.1303131   6.6925970   6.6271716  -4.8774107   6.8801270
                # (...)
                #           13
                #    1   1.7575832
                #    2   8.2313758
                #    3  -5.5451689
                elif line.startswith (" Calculating analytic gradient of the SCF energy"):
                    line = next (lines)
                    if line.startswith (" Gradient of AOints"):
                        nblocks = self.natoms / _ATOMS_PER_LINE
                        if (self.natoms % _ATOMS_PER_LINE != 0):
                            nblocks += 1
                        forces = []
                        for i in range (nblocks):
                            header      = next (lines)
                            coordinates = []
                            for j in range (3):
                                line   = next (lines)
                                tokens = TokenizeLine (line, converters=([None, ] + [float, ] * _ATOMS_PER_LINE))
                                coordinates.append (tokens[1:])
                            for (x, y, z) in zip (coordinates[0], coordinates[1], coordinates[2]):
                                if x != None:
                                    force = Force (
                                        x   =   -x   ,
                                        y   =   -y   ,
                                        z   =   -z   , )
                                    forces.append (force)
                        self.forces = forces
        except StopIteration:
            pass
        # . Close file
        lines.close ()


    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
