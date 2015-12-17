#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# . QM files
from MopacOutputFile     import MopacOutputFile
from GaussianOutputFile  import GaussianOutputFile
# . Molaris file
from MolarisAtomsFile    import MolarisAtomsFile

import subprocess, os.path


class QMCaller (object):
    """Base class to provide communication between Molaris and QM programs."""

    # . General options
    defaultAttributes = {
        "charge"               :   0              ,
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
        "cosmo"                :   True           ,
        # . If qmmm is True, fileAtoms should be set to "mol.in"
        "qmmm"                 :   False          ,
        "fileMopacError"       :   "run.err"      ,
        "fileMopacInput"       :   "run.mop"      ,
        "fileMopacOutput"      :   "run.out"      ,
        "pathMopac"            :   os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)

    def __init__ (self, **keywordArguments):
        super (QMCallerMopac, self).__init__ (**keywordArguments)


    def Run (self):
        # . Prepare a MOPAC calculation
        self.molaris.WriteMopacInput (filename=self.fileMopacInput, method=self.method, charge=self.charge, cosmo=self.cosmo)

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
        "fileGaussianError"    :   "job.err"      ,
        "fileGaussianInput"    :   "job.inp"      ,
        "fileGaussianOutput"   :   "job.log"      ,
        "pathGaussian"         :   os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)

    def __init__ (self, **keywordArguments):
        super (QMCallerGaussian, self).__init__ (**keywordArguments)


    def Run (self):
        # . Prepare a Gaussian calculation
        self.molaris.WriteGaussianInput (filename=self.fileGaussianInput, methodBasis=self.method, charge=self.charge)

        # . Run the calculation
        fileError  = open (self.fileGaussianError, "w")
        subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError)
        fileError.close ()

        # . Convert the output file from Gaussian to forces.out
        gaussian = GaussianOutputFile (filename=self.fileGaussianOutput)
        gaussian.WriteMolarisForces (filename=self.fileForces)
