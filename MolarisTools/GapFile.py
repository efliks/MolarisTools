#-------------------------------------------------------------------------------
# . File      : GapFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Utilities import TokenizeLine
import  collections


GapStep = collections.namedtuple ("GapStep", "reference  target")

_MODULE_LABEL = "GapFile"


class GapFile (object):
    """A class to represent a gap file from an evb_to_qm_map type of simulation."""

    def __init__ (self, filename, logging=True, ignoreStepZero=True):
        """Constructor."""
        self.logging        = logging
        self.ignoreStepZero = ignoreStepZero
        self._Parse (filename)


    @property
    def nsteps (self):
        if hasattr (self, "steps"):
            return len (self.steps)
        return 0


    def Extend (self, filename):
        """Extend with data from another file."""
        self._Parse (filename)


    def CalculateLRATerm (self, skip=None):
        """Calculate an LRA term, eg. <Eqmmm - Eevb>qmmm  or  <Eqmmm - Eevb>evb.

        Skip can be either the number (absolute) or the fraction (relative) of initial configurations to exclude from the calculation."""
        if   isinstance (skip, int  ):
            if self.logging:
                print ("# . %s> Skipping first %d configurations" % (_MODULE_LABEL, skip))
            steps = self.steps[skip:]
        elif isinstance (skip, float):
            nskip = int (self.nsteps * skip)
            if self.logging:
                percent = int (skip * 100.)
                print ("# . %s> Skipping first %d%% (%d) configurations" % (_MODULE_LABEL, percent, nskip))
            steps = self.steps[nskip:]
        else:
            steps = self.steps
        collect = 0.
        for (n, step) in enumerate (steps):
            collect += (step.target - step.reference)
        return (collect / n)


    def _Parse (self, filename):
        lines  = open (filename)
        steps  = []
        if self.logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, filename))
        try:
            #                   0  0.00100000 0.00    2 0.0000 0.0000   0   0.0000000  -1.0000000
            #       0       -7215.96    -1601300.60
            #       1       -7230.30    -1601303.26
            #   (...)
            # . Read and ignore header
            line   = next (lines)
            # . Read steps
            while True:
                line   = next (lines)
                tokens = TokenizeLine (line, converters=[int, float, float])
                index, potentialReference, potentialTarget = tokens
                step   = GapStep (
                    target    = potentialTarget     ,
                    reference = potentialReference  , )
                if index < 1:
                    if self.ignoreStepZero:
                        continue
                steps.append (step)
        except StopIteration:
            pass
        # . Finalize
        lines.close ()
        if hasattr (self, "steps"):
            self.steps.extend (steps)
        else:
            self.steps = steps
        if self.logging:
            nsteps = len (steps)
            print ("# . %s> Read %d steps" % (_MODULE_LABEL, nsteps))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
