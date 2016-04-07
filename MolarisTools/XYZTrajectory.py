#-------------------------------------------------------------------------------
# . File      : XYZTrajectory.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Atom      import Atom
from    Utilities import TokenizeLine, WriteData
import  math, exceptions


class XYZStep (object):
    """A class to handle a single step in an XYZ trajectory."""

    def __init__ (self, atoms, comment=""):
        """Constructor."""
        self.atoms   = atoms
        self.comment = comment


    @property
    def natoms (self):
        return len (self.atoms)


    def Write (self, filename="step.xyz", append=False):
        """Write a step."""
        header  = "%d\n%s\n" % (self.natoms, self.comment)
        data    = [header, ]
        for atom in self.atoms:
            data.append (atom.line)
        WriteData (data, filename=filename, append=append)


#===============================================================================
class XYZTrajectory (object):
    """A class to handle trajectories in the XYZ format (with varying number of atoms)."""

    def __init__ (self, filename="qm.xyz"):
        """Constructor."""
        self.filename = filename
        self._Parse ()


    def __getitem__ (self, index):
        """Return a step."""
        if abs (index) >= self.nsteps:
            raise exceptions.StandardError ("Index %d is out of range." % index)
        if index < 0:
            index = self.nsteps + index
        return self.steps[index]


    @property
    def nsteps (self):
        if hasattr (self, "steps"):
            return len (self.steps)
        else:
            return 0


    def _Parse (self):
        """Parse an XYZ file."""
        lines = open (self.filename)
        steps = []
        try:
            while True:
                line    = lines.next ()
                # . Read number of atoms
                natoms  = TokenizeLine (line, converters=[int])[0]
                # . Read comment
                line    = lines.next ()
                comment = line[:-1]
                # . Read atoms
                atoms   = []
                for i in range (natoms):
                    line   = lines.next ()
                    # . Detect format first
                    tokens = TokenizeLine (line)
                    if len (tokens) < 5:
                        # . Simple format
                        tokens = TokenizeLine (line, converters=[None, float, float, float])
                        atom   = Atom (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3])
                    else:
                        # . Extended format with forces and charges
                        tokens = TokenizeLine (line, converters=[None, float, float, float, float, float, float, float, float])
                        atom   = Atom (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], fx=tokens[4], fy=tokens[5], fz=tokens[6], fm=tokens[7], charge=tokens[8])
                    atoms.append (atom)
                # . Create a step and add it to the list of steps
                step = XYZStep (atoms=atoms, comment=comment)
                steps.append (step)
        except StopIteration:
            pass
        lines.close ()
        # . Finish up
        self.steps = steps


    def WriteAll (self, filename="traj.xyz", start=0, stop=-1):
        """Write all steps."""
        # . Reset the file
        fp = open (filename, "w")
        fp.close ()
        # . Write steps
        for step in self.steps[start:stop]:
            step.Write (filename=filename, append=True)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
