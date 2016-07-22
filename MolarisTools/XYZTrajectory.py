#-------------------------------------------------------------------------------
# . File      : XYZTrajectory.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine, WriteData
import  math, exceptions, collections


_FORMAT_SIMPLE    = "%2s   %8.3f   %8.3f   %8.3f\n"
_FORMAT_EXTENDED  = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n"


TrajStep         = collections.namedtuple ("TrajStep"  , "atoms  comment")

TrajAtom         = collections.namedtuple ("TrajAtom"  , "label  x  y  z")
TrajAtomExtended = collections.namedtuple ("TrajAtomExtended"  , "label  x  y  z  fx  fy  fz  fm  charge")

_MAX_ATOMS       = 1000
_MAX_STEPS       = 1000000


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
                        tokens = TokenizeLine (line, converters=[None, ] + [float, ] * 3)
                        atom   = TrajAtom (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3])
                    else:
                        # . Extended format with forces and charges
                        tokens = TokenizeLine (line, converters=[None, ] + [float, ] * 8)
                        atom   = TrajAtomExtended (label=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], fx=tokens[4], fy=tokens[5], fz=tokens[6], fm=tokens[7], charge=tokens[8])
                    atoms.append (atom)
                # . Create a step and add it to the list of steps
                step = TrajStep (atoms=atoms, comment=comment)
                steps.append (step)
        except StopIteration:
            pass
        lines.close ()
        # . Finish up
        self.steps = steps


    def Write (self, filename="traj.xyz", rangeSteps=(0, _MAX_STEPS)):
        """Write all steps."""
        openfile = open (filename, "w")
        # . Write steps
        start, stop = rangeSteps
        for step in self.steps[start:stop]:
            natoms  = len (step.atoms)
            openfile.write ("%d\n%s\n" % (natoms, step.comment))
            for atom in step.atoms:
                if   isinstance (atom, TrajAtom):
                    openfile.write (_FORMAT_SIMPLE   % atom.label, atom.x, atom.y, atom.z)
                elif isinstance (atom, TrajAtomExtended):
                    openfile.write (_FORMAT_EXTENDED % atom.label, atom.x, atom.y, atom.z, atom.fx, atom.fy, atom.fz, atom.fm, atom.charge)
        # . Close the file
        openfile.close ()


    def _GetHeader (self, rangeAtoms):
        line        = "#    "
        step        = self.steps[0]
        start, stop = rangeAtoms
        atoms       = step.atoms[start:stop]
        # . Replace atom serials by labels
        convert     = {}
        for iatom, atom in enumerate (atoms):
            line = "%s%9s" % (line, atom.label.center (9))
            convert[iatom + 1] = atom.label
        return (convert, line)


    def _WriteAtomicProperty (self, atomicProperty, filename, rangeAtoms, rangeSteps):
        convert, header = self._GetHeader (rangeAtoms)
        output      = open (filename, "w")
        output.write ("%s\n" % header)
        start, stop = rangeSteps
        steps       = self.steps[start:stop]
        start, stop = rangeAtoms
        for istep, step in enumerate (steps, 1):
            line    = "%4d" % istep
            for iatom, atom in enumerate (step.atoms[start:stop]):
                if isinstance (atom, TrajAtomExtended):
                    if   atomicProperty == "charge":
                        line     = "%s  %7.3f" % (line, atom.charge)
                    elif atomicProperty == "force" :
                        line     = "%s  %7.3f" % (line, atom.fm)
                    else:
                        raise exceptions.StandardError ("Unknown atomic property: %s" % atomicProperty)
                else:
                    raise exceptions.StandardError ("Atom has no properties.")
            output.write ("%s\n" % line)
        output.close ()
        return convert


    def WriteGnuplotCharges (self, filename="charges.dat", rangeAtoms=(0, _MAX_ATOMS), rangeSteps=(0, _MAX_STEPS)):
        """Write charges in a format suitable for viewing in Gnuplot."""
        return self._WriteAtomicProperty ("charge", filename, rangeAtoms, rangeSteps)


    def WriteGnuplotForces (self, filename="forces.dat", rangeAtoms=(0, _MAX_ATOMS), rangeSteps=(0, _MAX_STEPS)):
        """Write forces in a format suitable for viewing in Gnuplot."""
        return self._WriteAtomicProperty ("force", filename, rangeAtoms, rangeSteps)


    def BinCharges (self, rangeAtoms=(0, _MAX_ATOMS), rangeSteps=(0, _MAX_STEPS), sampling=0.1, limits=None):
        """For each atom, calculate the distribution of its charge."""
        start, stop = rangeSteps
        steps       = self.steps[start:stop]
        nsteps      = len (steps)
        step        = steps[0]
        start, stop = rangeAtoms
        natoms      = len (step.atoms[start:stop])
        # . For each atom, collect its charges throughout the trajectory
        atoms       = []
        for i in range (natoms):
            charges = []
            for step in steps:
                atom = step.atoms[start + i]
                if isinstance (atom, TrajAtomExtended):
                    charges.append (atom.charge)
                else:
                    raise exceptions.StandardError ("Atom has no charge property.")
            atoms.append (charges)
        # . Find extreme charges
        if limits:
            (minc, maxc) = limits
        else:
            minc =  999.
            maxc = -999.
            for atom in atoms:
                charges = atom
                mint    = min (charges)
                maxt    = max (charges)
                if mint < minc:
                    minc = mint
                if maxt > maxc:
                    maxc = maxt
        # . Calculate the number of bins and recalculate the sampling parameter
        spread   = maxc - minc
        nbins    = int (spread / sampling)
        sampling = spread / nbins
        # . For each bin, calculate its boundaries
        boundaries = []
        for i in range (nbins):
            left     = i * sampling + minc
            right    = left + sampling
            boundary = (left, right)
            boundaries.append (boundary)
        # . For each atom, calculate counts in each bin
        atomData = []
        for atom in atoms:
            charges = atom
            bins    = [0] * nbins
            for charge in charges:
                for iboundary, (left, right) in enumerate (boundaries):
                    if   iboundary < 1:
                        if charge >= left and charge <  right:
                            break
                    elif iboundary > (nbins - 2):
                        if charge >  left and charge <= right:
                            break
                    else:
                        if charge >  left and charge <  right:
                            break
                bins[iboundary] += 1
            atomData.append (bins)
        # . Finish up
        self.binsInfo = (natoms, minc, maxc, spread, nbins, sampling)
        self.bins     = (boundaries, atomData)


    def BinsWrite (self, filename="histogram.dat"):
        """Write histogram data to a file in a format suitable for Gnuplot."""
        if hasattr (self, "bins"):
            (boundaries, atomData) = self.bins
            (natoms, minc, maxc, spread, nbins, sampling) = self.binsInfo
        else:
            raise exceptions.StandardError ("First calculate bins.")
        # . Write header
        message = "natoms=%d, minc=%.3f, maxc=%.3f, spread=%.3f, nbins=%d, sampling=%f" % (natoms, minc, maxc, spread, nbins, sampling)
        lines   = ["# %s" % message, ]
        line    = "#" + 9 * " "
        for (serial, atom) in enumerate (self.steps[0].atoms, 1):
            label = "%s%d" % (atom.label, serial)
            line  = "%s %5s" % (line, label.center (5))
        lines.append (line)
        # . Write histogram
        for iboundary, (left, right) in enumerate (boundaries):
            charge = (left + right) / 2.
            line   = "%7.3f  " % charge
            for atom in atomData:
                counts = atom
                count  = counts[iboundary]
                line   = "%s  %4d" % (line, count)
            lines.append (line)
        # . Write file
        fo = open (filename, "w")
        for line in lines:
            fo.write (line + "\n")
        fo.close ()


    def BinsAssign (self):
        """For each atom, get a charge from the most populous bin."""
        if hasattr (self, "bins"):
            (boundaries, atomData) = self.bins
            (natoms, minc, maxc, spread, nbins, sampling) = self.binsInfo
        else:
            raise exceptions.StandardError ("First calculate bins.")
        # . Iterate atoms
        for atomSerial, (ad, atom) in enumerate (zip (atomData, self.steps[0].atoms), 1):
            maxCount = 0
            for (i, count) in enumerate (ad):
                if count > maxCount:
                    maxCount = count
                    maxi     = i
            (left, right) = boundaries[maxi]
            charge = (left + right) / 2.
            print ("%3d  %4s    %5.2f" % (atomSerial, atom.label, charge))


    def AverageCharges (self):
        """For each atom, calculate its average charge from the trajectory."""
        nsteps = len (self.steps)
        if nsteps > 0:
            natoms   = len (self.steps[0].atoms)
            averages = []
            for i in range (natoms):
                collect = []
                for step in self.steps:
                    atom = step.atoms[i]
                    collect.append (atom.charge)
                averages.append (sum (collect) / nsteps)
            self.averageCharges = averages


    def AverageChargesWrite (self, filename=""):
        """Write average charges."""
        if hasattr (self, "averageCharges"):
            toFile = False
            if filename != "":
                output = open (filename, "w")
                toFile = True
            step   = self.steps[0]
            for atomSerial, (atom, averageCharge) in enumerate (zip (step.atoms, self.averageCharges), 1):
                line = "%3d  %4s    %5.2f" % (atomSerial, atom.label, averageCharge)
                if toFile:
                    output.write (line + "\n")
                else:
                    print (line)
            if toFile:
                output.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
