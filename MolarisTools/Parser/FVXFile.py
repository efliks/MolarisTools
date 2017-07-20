#-------------------------------------------------------------------------------
# . File      : FVXFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections, math

from MolarisTools.Utilities  import TokenizeLine


FVXAtom = collections.namedtuple ("Atom"  , "serial  charge  fx  fy  fz  vx  vy  vz  x  y  z  fm  vm")

class FVXFile (object):
    """A class to represent data from a fvx.dat file."""

    def __init__ (self, filename="fvx.dat"):
        """Constructor."""
        self.filename = filename
        self._Parse ()


    @property
    def nsteps (self):
        return len (self.steps)


    def _Parse (self):
        data  = open (self.filename)
        steps = []
        step  = []
        try:
            while True:
                line   = data.next ()
                if line.count ("trajec"):
                    # . Found a new step
                    tokens   = TokenizeLine (line, converters=[int, None])
                    trajStep = tokens[0]
                    if step:
                        # . Save the previous step
                        steps.append (step)
                        step = []
                elif line.count ("atom"):
                    tokens     = TokenizeLine (line, converters=[None, None, int, float])
                    serial, charge = tokens[2], tokens[3]
                    # . Read forces
                    line       = data.next ()
                    tokens     = TokenizeLine (line, converters=[None, float, float, float])
                    fx, fy, fz = tokens[1:]
                    # . Read velocities
                    line       = data.next ()
                    tokens     = TokenizeLine (line, converters=[None, float, float, float])
                    vx, vy, vz = tokens[1:]
                    # . Read coordinates
                    line       = data.next ()
                    tokens     = TokenizeLine (line, converters=[None, float, float, float])
                    x, y, z    = tokens[1:]
                    # . Calculate the force    magnitude
                    fm = math.sqrt (fx * fx + fy * fy + fz * fz)
                    # . Calculate the velocity magnitude
                    vm = math.sqrt (vx * vx + vy * vy + vz * vz)
                    # . Create an atom and add it to the current step
                    atom = FVXAtom (serial=serial, charge=charge, fx=fx, fy=fy, fz=fz, vx=vx, vy=vy, vz=vz, x=x, y=y, z=z, fm=fm, vm=vm)
                    step.append (atom)
        except StopIteration:
            pass
        # . Close the file
        data.close ()
        # . Are there any steps left?
        if step:
            steps.append (step)
        self.steps = steps


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
