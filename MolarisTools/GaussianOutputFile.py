#-------------------------------------------------------------------------------
# . File      : GaussianOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine, WriteData
import collections, exceptions

Atom     = collections.namedtuple ("Atom"     , "symbol x y z charge")
Force    = collections.namedtuple ("Force"    , "x y z")
ScanStep = collections.namedtuple ("ScanStep" , "Efinal atoms forces charges espcharges")


class GaussianOutputFile (object):
    """A class to read a Gaussian output file."""

    def __init__ (self, filename="run_gauss.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        scan  = []
        try:
            while True:
                line = next (lines)
                # . Get the number of atoms and their coordinates
                if line.count ("Input orientation:") or line.count ("Z-Matrix orientation:"):
                    for skip in range (4):
                        next (lines)
                    natoms = 0
                    atoms  = []
                    while True:
                        line = next (lines)
                        if line.count ("----"):
                            break
                        tokens = TokenizeLine (line, converters=[int, int, int, float, float, float])
                        atomicNumber, x, y, z = tokens[1], tokens[3], tokens[4], tokens[5]
                        atom = Atom (symbol=atomicNumberToSymbol[atomicNumber], x=x, y=y, z=z, charge=0.)
                        atoms.append (atom)
                        natoms += 1
                    self.atoms = atoms

                # . Get the final energy (for semiempirical calculations)
                elif line.count ("NIter="):
                    tokens      = TokenizeLine (line, converters=[None, float, None, None])
                    self.Efinal = tokens[1] * HARTREE_TO_KCAL_MOL

                # . Get the final energy (for ab initio/DFT calculations)
                elif line.count ("SCF Done"):
                    tokens      = TokenizeLine (line, converters=[None, None, None, None, float, None, None, int, None])
                    self.Efinal = tokens[4] * HARTREE_TO_KCAL_MOL

                # . Get ESP charges
                elif line.count ("Charges from ESP fit"):
                    next (lines)
                    next (lines)
                    self.espcharges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.espcharges.append (charge)

                # . Get Mulliken charges
                elif line.count ("Mulliken atomic charges:"):
                    next (lines)
                    self.charges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.charges.append (charge)

                # . Get forces
                # . http://www.gaussian.com/g_tech/g_ur/k_force.htm
                # elif line.count ("***** Axes restored to original set *****"):
                #     for skip in range (4):
                #         next (lines)
                elif line.count ("Forces (Hartrees/Bohr)"):
                    next (lines)
                    next (lines)
                    self.forces = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, int, float, float, float])
                        force  = Force (x=tokens[2] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, y=tokens[3] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, z=tokens[4] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM)
                        self.forces.append (force)

                # . Get timing information
                elif line.count ("Job cpu time"):
                    tokens = TokenizeLine (line, converters=[None, None, None, int, None, int, None, int, None, float, None])
                    days, hours, minutes, seconds = tokens[3], tokens[5], tokens[7], tokens[9]
                    self.timings = {"days" : days, "hours" : hours, "minutes" : minutes, "seconds" : seconds}
                    # . Quit here, since composite jobs are not supported (?)
                    # break

                # . Determine if we have reached the end of an IRC step
                elif line.count ("-- Optimized point #"):
                    newStep = ScanStep (Efinal=self.Efinal, atoms=self.atoms[:], forces=self.forces[:], charges=self.charges[:], espcharges=[])  # <--FIX ME
                    scan.append (newStep)
        except StopIteration:
            pass
        # . Does the job involve a scan (IRC or PES)?
        if scan:
            self.scan = scan


    @property
    def nscan (self):
        if hasattr (self, "scan"):
            return len (self.scan)
        else:
            return 0

    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        else:
            return 0


    def WriteScanTrajectory (self, filename="traj.xyz", relative=True, reverse=False, append=False):
        if not hasattr (self, "scan"):
            raise exceptions.StandardError ("No trajectory found.")
        if relative:
            if reverse:
                Eref = self.scan[-1].Efinal
            else:
                Eref = self.scan[ 0].Efinal
        else:
            Eref = 0.
        if reverse:
            direction = -1
        else:
            direction =  1
        data = []
        for istep, step in enumerate (self.scan[::direction], 1):
            Erel = step.Efinal - Eref
            data.append ("%d\nStep %d: %f\n" % (self.natoms, istep, Erel))
            for atom in step.atoms:
                data.append ("%2s  %14.6f  %14.6f  %14.6f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        # . Write the file
        WriteData (data, filename, append=append)


    def WriteMolarisForces (self, filename="forces.out", Eref=0., useESPCharges=False):
        """Write a file in the Molaris-suitable format."""
        data = []
        data.append ("%f\n" % (self.Efinal - Eref))
        for force in self.forces:
            data.append ("%14.6f  %14.6f  %14.6f\n" % (force.x, force.y, force.z))
        # . Do not use CHELPG charges with semi-empirical methods?
        charges = self.charges
        if useESPCharges:
            charges = self.espcharges
        for charge in charges:
            data.append ("%14.6f\n" % charge)
        # . Write the file
        WriteData (data, filename)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
