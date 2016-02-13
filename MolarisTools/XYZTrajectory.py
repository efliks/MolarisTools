#-------------------------------------------------------------------------------
# . File      : XYZTrajectory.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine, WriteData
import  math, exceptions


class XYZAtom (object):
    """A class to handle an atom in an XYZ trajectory."""

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in keywordArguments.iteritems ():
            setattr (self, key, value)
        # label  x  y  z  fx  fy  fz  fm  charge

        # . Calculate the magnitude of the force
        checks = (
            not hasattr (self, "fm") ,
                hasattr (self, "fx") ,)
        if all (checks):
            self.fm = math.sqrt (self.fx ** 2 + self.fy ** 2 + self.fz ** 2)


#===============================================================================
_FORMAT_QM_EXT    = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.4f\n"
_FORMAT_QM_SIMPLE = "%2s   %8.3f   %8.3f   %8.3f\n"

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
            if hasattr (atom, "charge"):
                # . Write an atom with forces and a charge
                data.append (_FORMAT_QM_EXT % (atom.label, atom.x, atom.y, atom.z, atom.fx, atom.fy, atom.fz, atom.fm, atom.charge))
            else:
                # . Write a simple atom
                data.append (_FORMAT_QM_SIMPLE % (atom.label, atom.x, atom.y, atom.z))
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
                        atom   = XYZAtom (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3])
                    else:
                        # . Extended format with forces and charges
                        tokens = TokenizeLine (line, converters=[None, float, float, float, float, float, float, float, float])
                        atom   = XYZAtom (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], fx=tokens[4], fy=tokens[5], fz=tokens[6], fm=tokens[7], charge=tokens[8])
                    atoms.append (atom)
                # . Create a step and add it to the list of steps
                step = XYZStep (atoms=atoms, comment=comment)
                steps.append (step)
        except StopIteration:
            pass
        lines.close ()
        # . Finish up
        self.steps = steps


    def WriteAll (self, filename="traj.xyz"):
        """Write all steps."""
        # . Reset the file
        fp = open (filename, "w")
        fp.close ()
        # . Write steps
        for step in self.steps:
            step.Write (filename=filename, append=True)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
