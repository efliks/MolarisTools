#-------------------------------------------------------------------------------
# . File      : MolarisInputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine


class MolarisInputFile (object):
    """A class to represent a Molaris input file."""

    def __init__ (self, filename, pair):
        self.filename = filename
        self.pair     = pair
        self._Parse ()

    def _Parse (self):
        lines = open (self.filename)
        equil = None
        try:
            while True:
                line = lines.next ()
                # . Remove comments
                line = line[:line.find ("#")]

                if line.count ("constraint_pair"):
                    foo, aserial, bserial, forceConst, equilDist = TokenizeLine (line, converters=[None, int, int, float, float])
                    pair  = (aserial, bserial)
                    found = False
                    if pair == self.pair:
                        found = True
                    else:
                        pair = (bserial, aserial)
                        if pair == self.pair:
                            found = True
                    if found:
                        equil = [forceConst, equilDist]
                        break
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        self.equil = equil


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
