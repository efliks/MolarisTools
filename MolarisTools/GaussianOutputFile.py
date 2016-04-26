#-------------------------------------------------------------------------------
# . File      : GaussianOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine, WriteData
import collections, exceptions

Atom        = collections.namedtuple ("Atom"     , "symbol  x  y  z  charge")
Force       = collections.namedtuple ("Force"    , "x  y  z")
ScanStep    = collections.namedtuple ("ScanStep" , "Efinal  atoms  forces  charges  espcharges")

# . Coordinates, electric field, potential, charge
PointCharge = collections.namedtuple ("PointCharge" , "x  y  z  ex  ey  ez  potential  charge")


class GaussianOutputFile (object):
    """A class to read a Gaussian output file."""

    def __init__ (self, filename="run_gauss.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        scan  = []
        opt   = []
        try:
            while True:
                line = next (lines)
                # . Get the version and revision of Gaussian
                if line.startswith (" Gaussian"):
                    if line.count ("Revision"):
                        tokens = TokenizeLine (line, converters=[None, None, None, None])
                        self.version  = tokens[1][:-1]
                        self.revision = tokens[3][:-1]

                # . Get the number of atoms and their coordinates
                elif line.count ("Input orientation:") or line.count ("Z-Matrix orientation:"):
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
                # SCF Done:  E(RB+HF-LYP) =  -882.208703983     A.U. after   28 cycles
                elif line.count ("SCF Done"):
                    tokens      = TokenizeLine (line, converters=[None, None, None, None, float, None, None, int, None])
                    self.Efinal = tokens[4] * HARTREE_TO_KCAL_MOL


                # . Get the self energy of the charges
                # . In g03, there is no "a.u." at the end
                # Self energy of the charges =      -252.7809376522 a.u.
                elif line.count ("Self energy of the charges"):
                    tokens      = TokenizeLine (line, converters=[None] * 6 + [float])
                    self.Echrg  = tokens[-1] * HARTREE_TO_KCAL_MOL

                # . Get ESP charges (can be Merz-Kollman or CHELPG)
                elif line.count ("Charges from ESP fit"):
                    for i in range (2):
                        next (lines)
                    self.espcharges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.espcharges.append (charge)


                # . Get Mulliken charges
                # . The second condition is for Gaussian 09
                elif line.startswith (" Mulliken atomic charges:") or line.startswith (" Mulliken charges:"):
                    next (lines)
                    self.charges = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.charges.append (charge)


                # . Get forces
                # . http://www.gaussian.com/g_tech/g_ur/k_force.htm
                # . It looks like there is no need to multiply by -1, since Gaussian prints actual forces, not gradients (?)
                # elif line.count ("***** Axes restored to original set *****"):
                #     for skip in range (4):
                #         next (lines)
                elif line.count ("Center     Atomic                   Forces (Hartrees/Bohr)"):
                    for i in range (2):
                        next (lines)
                    self.forces = []
                    for i in range (natoms):
                        tokens = TokenizeLine (next (lines), converters=[int, int, float, float, float])
                        force  = Force (x=tokens[2] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, y=tokens[3] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, z=tokens[4] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM)
                        self.forces.append (force)


                # . Read coordinates and values of point charges
                # Point Charges:
                # XYZ=    2.0006    1.0001    0.0000 Q=    0.1110 A=    0.0000 R=    0.0000 C=    0.0000
                # XYZ=    2.0009    2.0911    0.0000 Q=   -0.3675 A=    0.0000 R=    0.0000 C=    0.0000
                # XYZ=    1.4863    2.4537    0.8897 Q=    0.1110 A=    0.0000 R=    0.0000 C=    0.0000
                #  (...)
                # Sum of input charges=            0.000000
                elif line.startswith (" Point Charges:"):
                    points = []
                    line   = next (lines)
                    while line.startswith (" XYZ="):
                        tokens   = TokenizeLine (line, converters=[None, float, float, float, None, float])
                        (x, y, z), charge = tokens[1:4], tokens[5]
                        point    = (x, y, z, charge)
                        points.append (point)
                        line     = next (lines)


                # . Get electrostatic properties
                # **********************************************************************
                #
                #            Electrostatic Properties Using The SCF Density
                #
                # **********************************************************************
                #
                #       Atomic Center    1 is at   3.665580  6.467202 12.974383
                #       Atomic Center    2 is at   4.909670  6.386763 13.616169
                #   (...)
                #      Read-in Center 2400 is at   5.504554 14.162232 26.811879
                #      Read-in Center 2401 is at   5.086579 15.682876 27.049785
                elif line.count ("Electrostatic Properties Using The SCF Density"):
                    pass
                #    pcPositions = []
                #    for i in range (3):
                #        next (lines)
                #    line = next (lines)
                #    while not line.count ("----"):
                #        x   =   float (line[32:42])
                #        y   =   float (line[42:52])
                #        z   =   float (line[52:62])
                #        position = (x, y, z)
                #        if line.startswith ("      Read-in Center"):
                #            # . Add a point charge
                #            pcPositions.append (position)
                #        else:
                #            # . Add a nucleus
                #            pass
                #        line = next (lines)


                #              Electrostatic Properties (Atomic Units)
                #
                # -----------------------------------------------------------------
                #    Center     Electric         -------- Electric Field --------
                #               Potential          X             Y             Z
                # -----------------------------------------------------------------
                #    1 Atom    -14.711204     -0.022648      0.000626     -0.009472
                #    2 Atom    -22.331530      0.084739      0.046163     -0.012921
                # (...)
                elif line.count ("Electrostatic Properties (Atomic Units)"):
                    pointsElectric = []
                    atomsElectric  = []
                    for i in range (5):
                        next (lines)
                    line = next (lines)
                    while not line.count ("----"):
                        tokens = TokenizeLine (line, converters=[float, float, float, float, None], reverse=True)
                        ez, ey, ex, potential = tokens[:4]
                        field  = (ex, ey, ez, potential)
                        if tokens[4] != "Atom":
                            # . Electrostatic potential and field on a point charge
                            pointsElectric.append (field)
                        else:
                            # . Electrostatic potential and field on a nucleus
                            atomsElectric.append (field)
                        line   = next (lines)
                    # . Save point charges
                    pointCharges = []
                    for (x, y, z, charge), (ex, ey, ez, potential) in zip (points, pointsElectric):
                        pc = PointCharge (
                            x   =  x   ,
                            y   =  y   ,
                            z   =  z   ,
                            # . Convert Eh/bohr to kcal/(mol*A)
                            ex  =  ex * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  ,
                            ey  =  ey * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  ,
                            ez  =  ez * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  ,
                            potential  =  potential * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM  ,
                            charge     =  charge    ,
                            )
                        pointCharges.append (pc)
                    self.pointCharges = pointCharges


                # . Get timing information
                elif line.count ("Job cpu time"):
                    tokens = TokenizeLine (line, converters=[None, None, None, int, None, int, None, int, None, float, None])
                    days, hours, minutes, seconds = tokens[3], tokens[5], tokens[7], tokens[9]
                    self.timings = {"days" : days, "hours" : hours, "minutes" : minutes, "seconds" : seconds}
                    # . Quit here, since composite jobs are not supported (?)
                    # break


                # . Check for a failed job
                elif line.count ("Error termination"):
                    raise exceptions.StandardError ("Error termination in Gaussian output file %s." % self.inputfile)

                # . Determine if we have reached the end of an IRC step
                elif line.count ("-- Optimized point #"):
                    newStep = ScanStep (Efinal=self.Efinal, atoms=self.atoms[:], forces=self.forces[:], charges=self.charges[:], espcharges=[])  # <--FIX ME
                    scan.append (newStep)

                # . Determine if we have reached the end of a geometry optimization step
                elif line.startswith (" Berny optimization."):
                    if hasattr (self, "Efinal"):
                        optStep = ScanStep (Efinal=self.Efinal, atoms=self.atoms[:], forces=self.forces[:], charges=self.charges[:], espcharges=[])  # <--FIX ME
                        opt.append (optStep)
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . Does the job involve a scan (IRC or PES)?
        if scan: self.scan = scan
        # . Does the job involve a geometry optimization?
        if opt:  self.opt  = opt


    @property
    def nscan (self):
        if hasattr (self, "scan"):
            return len (self.scan)
        return 0

    @property
    def nopt (self):
        if hasattr (self, "opt"):
            return len (self.opt)
        return 0

    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0

    @property
    def ncharges (self):
        if hasattr (self, "pointCharges"):
            return len (self.pointCharges)
        return 0


    # . Argument "trajectory" can be either a PES scan or a geometry optimization
    def _WriteTrajectory (self, trajectory, filename, relative=True, reverse=False, append=False):
        if relative:
            if reverse:
                Eref = trajectory[-1].Efinal
            else:
                Eref = trajectory[ 0].Efinal
        else:
            Eref = 0.
        if reverse:
            direction = -1
        else:
            direction =  1
        data = []
        for istep, step in enumerate (trajectory[::direction], 1):
            Erel = step.Efinal - Eref
            data.append ("%d\nStep %d: %f\n" % (self.natoms, istep, Erel))
            for atom in step.atoms:
                data.append ("%2s  %14.6f  %14.6f  %14.6f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        # . Write the file
        WriteData (data, filename, append=append)


    def WriteScanTrajectory (self, filename="traj.xyz", relative=True, reverse=False, append=False):
        if not hasattr (self, "scan"):
            raise exceptions.StandardError ("No trajectory found.")
        self._WriteTrajectory (self.scan, filename=filename, relative=relative, reverse=reverse, append=append)


    def WriteOptTrajectory (self, filename="opt_traj.xyz", relative=True, reverse=False, append=False):
        if not hasattr (self, "opt"):
            raise exceptions.StandardError ("No optimization trajectory found.")
        self._WriteTrajectory (self.opt, filename=filename, relative=relative, reverse=reverse, append=append)


    def WriteLastGeometry (self, filename="last_opt.xyz"):
        title = "%d\nLast step\n" % self.natoms
        data  = [title]
        for atom in self.atoms:
            data.append ("%2s  %14.6f  %14.6f  %14.6f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        # . Write the file
        WriteData (data, filename)


    def WriteMolarisForces (self, filename="forces.out", Eref=0., useESPCharges=False):
        """Write a file in the Molaris-suitable format."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
