#-------------------------------------------------------------------------------
# . File      : GAMESSOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine, WriteData
import collections, exceptions, datetime

Atom     = collections.namedtuple ("Atom"     , "symbol x y z charge")
Force    = collections.namedtuple ("Force"    , "x y z")


class GAMESSDatFile (object):
    """A class to read a GAMESS checkpoint file."""

    def __init__ (self, filename="run.dat"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                if line.startswith (" $VEC"):
                    vec = []
                    while True:
                        line = next (lines)
                        if line.startswith (" $END"):
                            break
                        else:
                            vec.append (line)
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . Save a list of lines (nothing more is needed)
        self.vec = vec


#-------------------------------------------------------------------------------
class GAMESSOutputFile (object):
    """A class to read a GAMESS output file."""

    def __init__ (self, filename="run.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get gradients
                #                         ----------------------
                #                         GRADIENT OF THE ENERGY
                #                         ----------------------
                #
                # UNITS ARE HARTREE/BOHR    E'X               E'Y               E'Z 
                #    1 P                0.000961748       0.034900593      -0.184607470
                #    2 O               -0.056832094       0.021746718      -0.036701933
                # ...
                if   line.count ("GRADIENT OF THE ENERGY"):
                    for i in range (3):
                        line = next (lines)
                    # . There are two sections that start from "GRADIENT OF THE ENERGY", we are interested in the second one
                    if line.count ("UNITS ARE HARTREE/BOHR"):
                        self.forces = []
                        while True:
                            tokens     = TokenizeLine (next (lines), converters=[int, None, float, float, float])
                            if not tokens:
                                break
                            fx, fy, fz = map (lambda convert: convert * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, tokens[2:5])
                            force      = Force (x=fx, y=fy, z=fz)
                            self.forces.append (force)


                # . Get Mulliken charges
                #          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
                #       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
                #    1 P            14.103960    0.896040        14.020001    0.979999
                #    2 O             8.406488   -0.406488         8.409966   -0.409966
                # ...
                elif line.count ("TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS"):
                    next (lines)
                    self.charges = []
                    while True:
                        tokens     = TokenizeLine (next (lines), converters=[int, None, float, float, float, float])
                        if not tokens:
                            break
                        mulliken   = tokens[3]
                        self.charges.append (mulliken)


                # . Get the final energy
                #                       TOTAL ENERGY =   -2179.8517336269
                elif line.count ("TOTAL ENERGY ="):
                    tokens      = TokenizeLine (line, converters=[None, None, None, float])
                    self.Efinal = tokens[3] * HARTREE_TO_KCAL_MOL


                # . Get atomic coordinates
                # ATOM      ATOMIC                      COORDINATES (BOHR)
                #           CHARGE         X                   Y                   Z
                # P          15.0     7.9538566823        2.4755410439       30.4000219645
                # O           8.0    10.6145908730        2.1127136543       31.2258322211
                # ...
                elif line.count ("COORDINATES (BOHR)"):
                    next (lines)
                    self.atoms = []
                    while True:
                        tokens = TokenizeLine (next (lines), converters=[None, float, float, float, float])
                        if not tokens:
                            break
                        atom   = Atom (symbol=tokens[0], x=tokens[2], y=tokens[3], z=tokens[4], charge=0.)
                        self.atoms.append (atom)


                # . Get timing information
                # TOTAL WALL CLOCK TIME=      102.7 SECONDS, CPU UTILIZATION IS  99.28%
                elif line.count ("TOTAL WALL CLOCK TIME"):
                    tokens       = TokenizeLine (line, converters=[None, None, None, None, float])
                    seconds      = tokens[-1]
                    hours, minutes, seconds = str (datetime.timedelta (seconds=seconds)).split (":")
                    self.timings = {"days" : 0, "hours" : hours, "minutes" : minutes, "seconds" : seconds} # . FIXME


                # . Get version number
                #          *         GAMESS VERSION =  5 DEC 2014 (R1)          *
                elif line.count ("GAMESS VERSION"):
                    tokens  = TokenizeLine (line, converters=[None] * 9)
                    version = " ".join ((tokens[4], tokens[5], tokens[6], tokens[7]))
                    self.version = version
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . TODO: Scans and optimizations as in Gaussian


    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        else:
            return 0


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
