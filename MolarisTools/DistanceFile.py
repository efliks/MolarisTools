#-------------------------------------------------------------------------------
# . File      : DistanceFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Utilities import TokenizeLine


class DistanceFile (object):
    """A class to represent data from a dist.dat file."""

    def __init__ (self, filename):
        self.filename = filename
        self.Parse_ToList ()
        self.Parse_ToDictionary ()

    @property
    def npairs (self):
        return len (self.pairs)

    @property
    def nsteps (self):
        return len (self.steps)


    def Parse_ToList (self):
        data  = open (self.filename)
        steps = []
        try:
            while True:
                line = data.next ()
                if line.count ("steps"):
                    pairs = []
                    while True:
                        line = data.next ()
                        if not (line.count ("average distance") or line.strip () == ""):
                            tokens = TokenizeLine (line, reverse=True, converters=[float, float, int, int, int])
                            electro, distance, pairb, paira, step = tokens
                            if step:
                                if pairs:
                                    steps.append (pairs)
                                    pairs = []
                            pairs.append (distance)
        except StopIteration:
            pass
        data.close ()
        if pairs:
            steps.append (pairs)
        self.steps  = steps


    def Parse_ToDictionary (self):
        data  = open (self.filename)
        pairs = {}
        try:
            while True:
                line = data.next ()
                if line.count ("steps"):
                    while True:
                        line = data.next ()
                        if not (line.count ("average distance") or line.strip () == ""):
                            tokens = TokenizeLine (line, reverse=True, converters=[float, float, int, int, int])
                            electro, distance, pairb, paira, step = tokens
                            key = (paira, pairb)
                            if not pairs.has_key (key):
                                pairs[key] = []
                            distances = pairs[key]
                            distances.append (distance)
        except StopIteration:
            pass
        data.close ()
        self.pairs = pairs


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
