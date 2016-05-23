#-------------------------------------------------------------------------------
# . File      : ORCAOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, HARTREE_TO_KCAL_MOL
from   Utilities import TokenizeLine
import collections, exceptions


Atom        = collections.namedtuple ("Atom"     , "symbol  x  y  z  charge")
Force       = collections.namedtuple ("Force"    , "x  y  z")

# . Convert from atomic units to Molaris units?
_DEFAULT_CONVERT_UNITS     = True

# . Multiply components of gradients by -1
_DEFAULT_REVERSE_GRADIENTS = False


class PCgradFile (object):
    """A class to read a file containing forces on point charges."""

    def __init__ (self, filename="run.pcgrad", reverse=_DEFAULT_REVERSE_GRADIENTS, convert=_DEFAULT_CONVERT_UNITS):
        """Constructor."""
        self.inputfile = filename
        self._Parse (reverse=reverse, convert=convert)

    def _Parse (self, reverse, convert):
        lines = open (self.inputfile)
        # . Get the number of atoms
        line   = next (lines)
        natoms = int (line)
        # . Read forces
        forces = []
        for i in range (natoms):
            line   = next (lines)
            gx, gy, gz = TokenizeLine (line, converters=[float, float, float])
            if reverse:
                gx, gy, gz = (-gx, -gy, -gz)
            if convert:
                gx, gy, gz = (gx * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, gy * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, gz * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM)
            force = Force (
                x   =   gx  ,
                y   =   gy  ,
                z   =   gz  ,
                )
            forces.append (force)
        self.forces = forces
        # . Close the file
        lines.close ()


#-------------------------------------------------------------------------------
class EngradFile (object):
    """A class to read a file containing forces on quantum atoms."""

    def __init__ (self, filename="run.engrad", reverse=_DEFAULT_REVERSE_GRADIENTS, convert=_DEFAULT_CONVERT_UNITS):
        """Constructor."""
        self.inputfile = filename
        self._Parse (reverse=reverse, convert=convert)

    def _Parse (self, reverse, convert):
        lines = open (self.inputfile)
        for i in range (3):
            next (lines)
        # . Get the number of atoms
        line   = next (lines)
        natoms = int (line)
        # . Skip a few lines
        for i in range (7):
            next (lines)
        # . Read forces
        forces = []
        for i in range (natoms):
            templ = []
            for j in range (3):
                line = next (lines)
                component = float (line)
                if reverse:
                    component = -component
                if convert:
                    component =  component * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM 
                templ.append (component)
            force = Force (
                x   =   templ[0]    ,
                y   =   templ[1]    ,
                z   =   templ[2]    ,
                )
            forces.append (force)
        self.forces = forces
        # . Close the file
        lines.close ()
    

#-------------------------------------------------------------------------------
class ORCAOutputFile (object):
    """A class to read an ORCA output file."""

    def __init__ (self, filename="run.out", reverse=_DEFAULT_REVERSE_GRADIENTS, convert=_DEFAULT_CONVERT_UNITS):
        """Constructor."""
        self.inputfile = filename
        self._Parse (reverse=reverse, convert=convert)


    def _Parse (self, reverse, convert):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get coordinates of QM atoms
                # ---------------------------------
                # CARTESIAN COORDINATES (ANGSTROEM)
                # ---------------------------------
                #   C      5.663910    4.221157   -1.234141
                #   H      5.808442    3.140412   -1.242145
                # (...)
                if line.startswith ("CARTESIAN COORDINATES (ANGSTROEM)"):
                    next (lines)
                    line     = next (lines)
                    geometry = []
                    while line != "\n":
                        tokens = TokenizeLine (line, converters=[None, float, float, float])
                        label, x, y, z = tokens
                        atom   = (label, x, y, z)
                        geometry.append (atom)
                        line   = next (lines)


                # . Get charges on QM atoms
                # -----------------------
                # MULLIKEN ATOMIC CHARGES
                # -----------------------
                #    0 C :   -0.520010
                #    1 H :    0.271953
                # (...)
                # Sum of atomic charges:   -0.0000000
                if line.startswith ("MULLIKEN ATOMIC CHARGES"):
                    next (lines)
                    line    = next (lines)
                    charges = []
                    while not line.startswith ("Sum of atomic charges:"):
                        tokens = TokenizeLine (line, separator=":", converters=[None, float])
                        charge = tokens[1]
                        charges.append (charge)
                        line   = next (lines)
                    self.charges = charges
                    # . Construct the final list of atoms with charges
                    atoms   = []
                    for (label, x, y, z), charge in zip (geometry, charges):
                        atom   = Atom (
                            symbol  =   label   ,
                            x       =   x       ,
                            y       =   y       ,
                            z       =   z       ,
                            charge  =   charge  ,
                            )
                        atoms.append (atom)
                    self.atoms = atoms


                # . Get gradients on QM atoms
                # ------------------
                # CARTESIAN GRADIENT
                # ------------------
                # 
                #    1   C   :   -0.017273415    0.000431161    0.011902545
                #    2   H   :    0.011246801   -0.004065387   -0.003492146
                # (...)
                elif line.startswith ("CARTESIAN GRADIENT"):
                    for i in range (2):
                        next (lines)
                    line   = next (lines)
                    forces = []
                    while line != "\n":
                        tokens = TokenizeLine (line, converters=[int, None, None, float, float, float])
                        gx, gy, gz = tokens[3:6]
                        if reverse:
                            gx, gy, gz = (-gx, -gy, -gz)
                        if convert:
                            gx, gy, gz = (gx * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, gy * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, gz * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM)
                        force  = Force (
                            x   =   gx  ,
                            y   =   gy  ,
                            z   =   gz  ,
                            )
                        forces.append (force)
                        line   = next (lines)
                    self.forces = forces


                # . Get the final energy
                # FINAL SINGLE POINT ENERGY      -263.834308915009
                elif line.startswith ("FINAL SINGLE POINT ENERGY"):
                    tokens = TokenizeLine (line, converters=[float, ], reverse=True)
                    if convert:
                        self.Efinal = tokens[-1] * HARTREE_TO_KCAL_MOL
                    else:
                        self.Efinal = tokens[-1]
        except StopIteration:
            pass
        # . Close the file
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
