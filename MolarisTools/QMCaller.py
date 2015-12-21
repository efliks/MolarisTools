#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# . QM files
from GaussianInputFile   import GaussianInputFile
from GaussianOutputFile  import GaussianOutputFile
from MopacOutputFile     import MopacOutputFile
# . Molaris file
from MolarisAtomsFile    import MolarisAtomsFile

import subprocess, os.path, exceptions


class QMCaller (object):
    """Base class to provide communication between Molaris and QM programs."""

    # . General options
    defaultAttributes = {
        "charge"               :   0              ,
        "multiplicity"         :   1              ,
        "method"               :   "PM3"          ,
        "fileAtoms"            :   "atoms.inp"    ,
        "fileForces"           :   "forces.out"   ,
        "fileTrajectory"       :   "qm.xyz"       ,
            }

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems ():                 setattr (self, key, value)
        # . Preparatory step
        self._Prepare ()


    def _Prepare (self):
        # . Read atoms.inp from Molaris
        self.molaris = MolarisAtomsFile (filename=self.fileAtoms)

        # . Write geometries for the purpose of viewing
        if self.fileTrajectory:
            self.molaris.WriteQM (filename=self.fileTrajectory, link=True, append=True)


#===============================================================================
class QMCallerMopac (QMCaller):
    """A class to provide communication between Molaris and Mopac."""

    # . Options specific to Mopac
    defaultAttributes = {
        "cosmo"                :   False          ,
        # . If qmmm is True, Mopac will use point charges to polarize the wavefunction
        # . fileAtoms should be set to "mol.in"
        "qmmm"                 :   False          ,
        "fileMopacError"       :   "run.err"      ,
        "fileMopacInput"       :   "run.mop"      ,
        "fileMopacOutput"      :   "run.out"      ,
        "pathMopac"            :   os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)

    def __init__ (self, **keywordArguments):
        super (QMCallerMopac, self).__init__ (**keywordArguments)

        if self.cosmo and self.qmmm:
            raise exceptions.StandardError ("Both cosmo and qmmm options cannot be enabled.")
        if self.qmmm:
            if self.fileAtoms != "mol.in":
                raise exceptions.StandardError ("With qmmm option on, fileAtoms can only be mol.in.")


    def Run (self):
        # . Prepare a MOPAC calculation
        self.molaris.WriteMopacInput (filename=self.fileMopacInput, method=self.method, charge=self.charge, cosmo=self.cosmo, qmmm=self.qmmm)

        # . Run the calculation
        fileError  = open (self.fileMopacError, "w")
        subprocess.check_call ([self.pathMopac, self.fileMopacInput], stdout=fileError, stderr=fileError)
        fileError.close ()

        # . Convert the output file from MOPAC to forces.out
        mopac = MopacOutputFile (filename=self.fileMopacOutput)
        mopac.WriteMolarisForces (filename=self.fileForces)


#===============================================================================
class QMCallerGaussian (QMCaller):
    """A class to provide communication between Molaris and Gaussian."""

    # . Options specific to Gaussian
    defaultAttributes = {
        # . env may define variables such as GAUSS_EXEDIR and GAUSS_SCRDIR
        "env"                     :   None         ,
        "ncpu"                    :   1            ,
        "memory"                  :   1            ,
        # . If qmmm is True, Gaussian will use point charges to polarize the wavefunction
        "qmmm"                    :   False        ,
        # . restart takes the wavefunction from the previous calculation
        "restart"                 :   False        ,
        "fileGaussianError"       :   "job.err"    ,
        "fileGaussianInput"       :   "job.inp"    ,
        "fileGaussianOutput"      :   "job.log"    ,
        "fileGaussianCheckpoint"  :   "job.chk"    ,
        "pathGaussian"            :   os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)

    def __init__ (self, **keywordArguments):
        super (QMCallerGaussian, self).__init__ (**keywordArguments)

        # . Determine if a semiempirical potential is used
        method = self.method[:3]
        if method in ("AM1", "PM3", ) and self.qmmm:
            raise exceptions.StandardError ("Point charges cannot be used with semiempirical methods.")


    def Run (self):
        # . Prepare a Gaussian calculation
        gaussian = GaussianInputFile (
            qmmm           =  self.qmmm                    ,
            ncpu           =  self.ncpu                    ,
            memory         =  self.memory                  ,
            method         =  self.method                  ,
            charge         =  self.charge                  ,
            multiplicity   =  self.multiplicity            ,
            fileInput      =  self.fileGaussianInput       ,
            fileCheckpoint =  self.fileGaussianCheckpoint  ,
            restart        =  os.path.exists (self.fileGaussianCheckpoint) and self.restart ,
                )
        atoms        = (self.molaris.qatoms + self.molaris.latoms)
        pointCharges = (self.molaris.patoms + self.molaris.watoms)
        gaussian.Write (atoms, pointCharges)

        # . Run the calculation
        fileError  = open (self.fileGaussianError, "w")
        if self.env:
            subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError, env=self.env)
        else:
            subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError)
        fileError.close ()

        # . Convert the output file from Gaussian to forces.out
        output = GaussianOutputFile (filename=self.fileGaussianOutput)
        output.WriteMolarisForces (filename=self.fileForces)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
