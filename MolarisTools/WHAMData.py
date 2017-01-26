#-------------------------------------------------------------------------------
# . File      : WHAMData.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Utilities import TokenizeLine
import glob

_DefaultRAMInput = \
"""! wham input file
%f        ! dz
%f %f     ! zmin zmax
0.001     ! tolw
100       ! testfrq
.false.   ! pmfmon
0         ! nskip
.false.   ! periodic
.false.   ! Use degree/kcal/mol/rad^2, NOT Angstrom,kcal/mol/A^2
.false.   ! Symmetrize biased density?
.false.   ! don't wrap angles?
-1
""" 


class WHAMTwoDistances (object):
    """A class to read Molaris files and prepare input files for the WHAM program."""

    # . Leave and attack define pairs of atoms for which reaction coordinate is calculated
    def __init__ (self, nwindows=0, leave=(1, 2), attack=(2, 6), nsteps=5000, timestep=0.001, qmmmInterval=10):
        if nwindows < 1:
            dirs     = glob.glob ("w[0-9]*")
            nwindows = len (dirs)
            print ("Found %d windows." % nwindows)
        self.leave    = leave
        self.attack   = attack
        self.nwindows = nwindows
        self.nsteps       = nsteps
        self.timestep     = timestep
        self.qmmmInterval = qmmmInterval
        self._Parse ()


    def _Parse (self):
        equil   = []
        windows = []
        for w in range (1, self.nwindows + 1):
            distFile     = DistanceFile ("w%d/window_%d/dist.dat" % (w, w))
            distLeave    = distFile.pairs[self.leave]
            distAttack   = distFile.pairs[self.attack]
            reactionCoor = []
            for dl, da in zip (distLeave, distAttack):
                rc = dl - da
                reactionCoor.append (rc)
            averageCoor = sum (reactionCoor) / len (reactionCoor)

            inputfile = "w%d/window_%d.inp" % (w, w)
            lines     = open (inputfile)
            try:
                while True:
                    line = lines.next ()
                    if   line.count ("constraint_pair"):
                        foo, aserial, bserial, forceConst, equilDist = TokenizeLine (line, converters=[None, int, int, float, float])
                        pair    = (aserial, bserial)
                        pairRev = (bserial, aserial)
                        if   (self.leave  == pair) or (self.leave  == pairRev):
                            forceConstLeave  = forceConst
                            equilDistLeave   = equilDist
                        elif (self.attack == pair) or (self.attack == pairRev):
                            forceConstAttack = forceConst
                            equilDistAttack  = equilDist
            except StopIteration:
                pass
            lines.close ()

            # reactionCoor, averageCoor, forceConstLeave, equilDistLeave, forceConstAttack, equilDistAttack
            equil.append ([forceConstLeave, averageCoor])
            windows.append (reactionCoor)
        self.equil     = equil
        self.windows   = windows


    def WriteMetaData (self, filename="wham.meta"):
        """Generate the metadata file."""
        metafile = open (filename, "w")
        # . DEBUG
        # print ("len (equil) = %d" % len (self.equil))
        for iwindow, (forceConst, equilDist) in enumerate (self.equil, 1):
            tsfile = "wham%02d.ts" % (iwindow)
            # . Multiply the force constant by 2 because Molaris used the constraining potential of a form:
            #       V = k (x-x0)**2, thus no 1/2 factor.
            metafile.write ("%s        %8.3f       %8.3f\n" % (tsfile, equilDist, forceConst * 2.))
        metafile.close ()


    def WriteTimeSeries (self):
        """Generate timeseries files."""
        for iwindow, window in enumerate (self.windows, 1):
            filename  = "wham%02d.ts" % iwindow
            output    = open (filename, "w")
            timeshift = (iwindow - 1) * self.timestep * self.qmmmInterval * self.nsteps
        
            # for idistance, (distance, Epot) in enumerate (zip (window, self.energies), 1):
                # time = timestep * (idistance - 1) * qmmmInterval + timeshift
                # time = timestep * idistance * qmmmInterval
                # output.write ("%4d      %8.3f       %8.3f       %8.3f\n" % (idistance, time, distance, Epot))
            for idistance, distance in enumerate (window, 1):
                output.write ("%4d      %8.3f\n" % (idistance, distance))
            output.close ()
            print ("Time-series file %s done." % filename)


#-------------------------------------------------------------------------------
class WHAMData (object):
    """A class to read Molaris files and prepare input files for the WHAM program."""

    def __init__ (self, rootDir=".", nwindows=21, nsteps=5000, timestep=0.001, qmmmInterval=10, reverse=False):
        self.rootDir      = rootDir
        self.nwindows     = nwindows
        self.nsteps       = nsteps
        self.timestep     = timestep
        self.qmmmInterval = qmmmInterval
        self._Parse ()
        if reverse:
            print ("Reversing the PMF ...")
            self.equil.reverse ()
            self.windows.reverse ()


    def _Parse (self):
        print ("Parsing input files ...")
        self._ReadInputs    ()
        print ("Parsing dist.dat files ...")
        self._ReadDistances ()
        print ("Parsing output files for potential energies ...")
        self._ReadPotentialEnergies ()


    def _ReadDistances (self):
        windows = []
        for w in range (1, self.nwindows + 1):
            distfile = "%s/w%d/window_%d/dist.dat" % (self.rootDir, w, w)
            lines    = iter (open (distfile).readlines ())
            try:
                while True:
                    line = lines.next ()
                    if line.count ("steps"):
                        distances = []
                        while True:
                            line = lines.next ()
                            if line.count ("average distance"):
                                break
                            tokens = TokenizeLine (line, reverse=True, converters=[float, float, int, int, int])
                            electro, distance, pairb, paira, step = tokens

                            distances.append (distance)
                        windows.append (distances)
            except StopIteration:
                pass
        self.windows = windows


    def _ReadInputs (self):
        equil = []
        for w in range (1, self.nwindows + 1):
            inputfile = "%s/w%d/window_%d.inp" % (self.rootDir, w, w)
            lines     = iter (open (inputfile).readlines ())
            tmpl      = []
            try:
                while True:
                    line = lines.next ()
                    if   line.count ("constraint_pair"):
                        entry, aserial, bserial, forceConst, equilDist = TokenizeLine (line, converters=[None, int, int, float, float])
                        tmpl.append ([aserial, bserial, forceConst, equilDist])

                    # . There can only be one such line because one cannot define a filename different than "dist.dat"
                    elif line.count ("dist_atoms"):
                        entry, adist, bdist = TokenizeLine (line, converters=[None, int, int])
            except StopIteration:
                pass
            # . Find the constraint_pair of interest
            for aserial, bserial, forceConst, equilDist in tmpl:
                if (aserial, bserial) == (adist, bdist) or (aserial, bserial) == (bdist, adist):
                    equil.append ([forceConst, equilDist])
                    break
        self.equil     = equil
        foo, distFirst = equil[ 0]
        foo, distLast  = equil[-1]
        binsize      = (distLast - distFirst) / self.nwindows
        self.binsize = binsize
        self.zmin    = distFirst - binsize / 2.
        self.zmax    = distLast  + binsize / 2.


    def _ReadPotentialEnergies (self):
        for iwindow, window in enumerate (self.windows, 1):
            # . Read Molaris output file to find Etot for each configuration
            molaris  = MolarisOutputFile (filename="%s/w%d/window_%d.out" % (self.rootDir, iwindow, iwindow))
            # . There is only one step per file, because here FEP is split into separate windows
            mold     = molaris.fepSteps[0]
            energies = []
            for moldStep in mold:
                if hasattr (moldStep, "classic"):
                    # . Take only energies calculated by QM-calls
                    energies.append (moldStep.system.epot)
        self.energies = energies


    def WriteMetaData (self, filename="wham.meta"):
        """Generate the metadata file."""
        metafile = open (filename, "w")
        # . DEBUG
        # print ("len (equil) = %d" % len (self.equil))
        for iwindow, (forceConst, equilDist) in enumerate (self.equil, 1):
            tsfile = "wham%02d.ts" % (iwindow)
            # . Multiply the force constant by 2 because Molaris used the constraining potential of a form:
            #       V = k (x-x0)**2, thus no 1/2 factor.
            metafile.write ("%s        %8.3f       %8.3f\n" % (tsfile, equilDist, forceConst * 2.))
        metafile.close ()


    def WriteTimeSeries (self):
        """Generate timeseries files."""
        for iwindow, window in enumerate (self.windows, 1):
            filename  = "wham%02d.ts" % iwindow
            output    = open (filename, "w")
            timeshift = (iwindow - 1) * self.timestep * self.qmmmInterval * self.nsteps
        
            for idistance, (distance, Epot) in enumerate (zip (window, self.energies), 1):
                # time = timestep * (idistance - 1) * qmmmInterval + timeshift
                # time = timestep * idistance * qmmmInterval
                # output.write ("%4d      %8.3f       %8.3f       %8.3f\n" % (idistance, time, distance, Epot))
                output.write ("%4d      %8.3f\n" % (idistance, distance))
            output.close ()
            print ("Time-series file %s done." % filename)


    def WriteTimeSeriesRAM (self, reverse=False, average=False):
        direction = 1 if not reverse else -1

        for iwindow, (window, (forceConst, equilDist)) in enumerate (zip (self.windows[::direction], self.equil[::direction]), 1):
            filename  = "pmf_%d_1.dat" % iwindow
            output    = open (filename, "w")
            # . If true, replace a declared distance by an average distance.
            if average:
                distance = sum (window) / len (window)
            else:
                distance = equilDist
            output.write ("# DROFF = %14.8f K = %14.8f\n" % (distance, forceConst * 2.))

            for distance in window:
                output.write ("%.3f\n" % distance)
            output.close ()
            print ("Wrote file %s." % filename)
#   # DROFF =     1.50000000 K =   200.00000000 
#   1.488
#   1.487


    def WriteInputRAM (self, filename="wham.inp", reverse=False):
        if reverse:
            binsize    = -self.binsize
            zmin, zmax =  self.zmax, self.zmin
        else:
            binsize    =  self.binsize
            zmin, zmax =  self.zmin, self.zmax
        output  = open (filename, "w")
        output.write (_DefaultRAMInput % (binsize, zmin, zmax))
        output.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
