#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : prep_wham_files.py
# . Copyright : USC, Mikolaj J. Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import exceptions, collections, math

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


Protein   =  collections.namedtuple ("Protein"  ,  " ebond    ethet     ephi    eitor    evdw     emumu     ehb_pp  ")
Water     =  collections.namedtuple ("Water"    ,  " ebond    ethet                      evdw     emumu     ehb_ww  ")
Prowat    =  collections.namedtuple ("Prowat"   ,  "                                     evdw     emumu     ehb_pw  ")
Long      =  collections.namedtuple ("Long"     ,  " elong                                                          ")
Ac        =  collections.namedtuple ("Ac"       ,  " evd_ac   evd_acw   ehb_ac  ehb_acw  emumuac  emumuacw          ")
Evb       =  collections.namedtuple ("Evb"      ,  " ebond    ethet     ephi             evdw     emumu     eoff  egaschift  eindq   ebulk ")
Induce    =  collections.namedtuple ("Induce"   ,  " eindp    eindw             ")
Const     =  collections.namedtuple ("Const"    ,  " ewatc    eproc     edistc  ")
Langevin  =  collections.namedtuple ("Langevin" ,  " elgvn    evdw_lgv  eborn   ")
Classic   =  collections.namedtuple ("Classic"  ,  " classic  quantum           ")
System    =  collections.namedtuple ("System"   ,  " epot     ekin      etot    ")

#-------------------------------------------------------------------------------
class MDStep (object):
    """A class to hold energies of an MD step."""

    def __init___ (self):
        # . Seems to have no effect anyway?
        for att in ("protein", "water", "prowat", "elong", "ac", "evb", "induce", "const", "langevin", "classic", "system", ):
            setattr (self, att, None)


#-------------------------------------------------------------------------------
class MolarisOutputFile (object):
    """A class for reading output files from Molaris."""

    def __init__ (self, filename="rs_fep.out"):
        """Constructor."""
        self.filename      = filename
        self.currentMDStep = None
        # . fepSteps are the FEP steps (usually 11), each consists of many MD steps (usually 500)
        self.fepSteps      = []
        self._Parse ()


    # . Returns the number of FEP steps (lambda 0 ... 1)
    @property
    def nfepSteps (self):
        return len (self.fepSteps)

    # . Returns the total number of MD steps
    @property
    def nmdSteps (self):
        return sum (map (lambda fepStep: len (fepStep), self.fepSteps))


    def _ReadProtein (self, line, lines):
        toka = line.split ()
        tokb = lines.next ().split ()
        tokc = lines.next ().split ()
        tokd = lines.next ().split ()
        protein = Protein (
                ebond  = float ( toka[4] ) ,
                ethet  = float ( toka[7] ) ,
                ephi   = float ( tokb[2] ) ,
                eitor  = float ( tokb[5] ) ,
                evdw   = float ( tokc[2] ) ,
                emumu  = float ( tokc[5] ) ,
                ehb_pp = float ( tokd[2] ) ,
                          )
        self.currentMDStep.protein = protein


    def _ReadSystem (self, line, lines):
        toka   = line.split ()
        system = System (
                epot = float ( toka[4]  ) ,
                ekin = float ( toka[7]  ) ,
                etot = float ( toka[10] ) ,
                        )
        self.currentMDStep.system = system


    def _ReadClassic (self, line, lines):
        toka     = line.split ()
        energies = Classic (
                classic = float ( toka[4] ) ,
                quantum = float ( toka[7] ) ,
                           )
        self.currentMDStep.classic = energies


    def _Parse (self):
        lines   = iter (open (self.filename).readlines ())
        mdSteps = []
        try:
            while True:
                line = lines.next ()
                if line.startswith (" Energies for the system at step"):
                    self.currentMDStep = MDStep ()
                    while True:
                        line = lines.next ()
                        if   line.startswith ( " protein"  ):
                            self._ReadProtein (line, lines)
                        elif line.startswith ( " water"    ):
                            pass
                        elif line.startswith ( " pro-wat"  ):
                            pass
                        elif line.startswith ( " long"     ):
                            pass
                        elif line.startswith ( " ac"       ):
                            pass
                        elif line.startswith ( " evb"      ):
                            pass
                        elif line.startswith ( " induce"   ):
                            pass
                        elif line.startswith ( " const."   ):
                            pass
                        elif line.startswith ( " langevin" ):
                            pass
                        elif line.startswith ( " classic"  ):
                            self._ReadClassic (line, lines)
                        elif line.startswith ( " system"   ):
                            self._ReadSystem (line, lines)
                            break
                    mdSteps.append (self.currentMDStep)

                elif line.startswith (" Average energies for the system at the step"):
                    self.fepSteps.append (mdSteps)
                    mdSteps = []
        except StopIteration:
            pass


#-------------------------------------------------------------------------------
def TokenizeLine (line, converters=None, separator=None, reverse=False):
    """Tokenize a line with optional converters and separator."""
    tokens = None
    if line is not None:
        if separator is None:
            tokens = line.split ()
        else:
            tokens = line.split (separator)
        # . The last token becomes first
        if reverse:
            tokens.reverse ()
        if converters is not None:
            # . Truncate tokens if there is fewer of them than converters
            ntokens     = len (tokens)
            nconverters = len (converters)
            if ntokens > nconverters:
                tokens = tokens[:nconverters]
            for (i, (token, converter)) in enumerate (zip (tokens, converters)):
                if converter is None:
                    new = token
                else:
                    try:
                        new = converter (token)
                    except:
                        raise exceptions.StandardError ("Unable to convert token " + repr (i) + ".", True)
                tokens[i] = new
    return tokens


#-------------------------------------------------------------------------------
class WHAMData (object):
    """A class to read Molaris files and prepare input files for the WHAM program."""

    def __init__ (self, nwindows=21, nsteps=5000, timestep=0.001, qmmmInterval=10):
        self.nwindows     = nwindows
        self.nsteps       = nsteps
        self.timestep     = timestep
        self.qmmmInterval = qmmmInterval
        self._Parse ()


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
            distfile = "w%d/window_%d/dist.dat" % (w, w)
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
                            tokens   = line.split ()
                            distance = float (tokens[3])
                            distances.append (distance)
                        windows.append (distances)
            except StopIteration:
                pass
        self.windows = windows


    def _ReadInputs (self):
        equil = []
        for w in range (1, self.nwindows + 1):
            inputfile = "w%d/window_%d.inp" % (w, w)
            lines     = iter (open (inputfile).readlines ())
            try:
                while True:
                    line = lines.next ()
                    if line.count ("constraint_pair"):
                        entry, aserial, bserial, forceConst, equilDist = TokenizeLine (line, converters=[None, int, int, float, float])
                        equil.append ([forceConst, equilDist])
            except StopIteration:
                pass
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
            molaris  = MolarisOutputFile (filename="w%d/window_%d.out" % (iwindow, iwindow))
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


    def WriteTimeSeriesRAM (self, reverse=False):
        direction = 1 if not reverse else -1

        for iwindow, (window, (forceConst, equilDist)) in enumerate (zip (self.windows[::direction], self.equil[::direction]), 1):
            filename  = "pmf_%d_1.dat" % iwindow
            output    = open (filename, "w")
            output.write ("# DROFF = %14.8f K = %14.8f\n" % (equilDist, forceConst * 2.))

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
wham = WHAMData (nwindows=31)

# . Write files in the original format of WHAM
wham.WriteMetaData   ()
wham.WriteTimeSeries ()

# . Write files in Ram's format
wham.WriteInputRAM      (reverse=True)
wham.WriteTimeSeriesRAM (reverse=True)
