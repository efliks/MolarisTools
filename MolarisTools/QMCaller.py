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

import exceptions, subprocess, os.path


class QMCaller (object):
    """A class to provide communication between Molaris and a QM program (Mopac or Gaussian)."""

    defaultAttributes = {
        "cosmo"                :   True           ,
        "charge"               :   0              ,
        "method"               :   "PM3"          ,
        "program"              :   "mopac"        ,
        "fileAtoms"            :   "atoms.inp"    ,
        "fileForces"           :   "forces.out"   ,
        "fileTrajectory"       :   "qm.xyz"       ,
        "fileMopacError"       :   "run.err"      ,
        "fileMopacInput"       :   "run.mop"      ,
        "fileMopacOutput"      :   "run.out"      ,
        "fileGaussianError"    :   "job.err"      ,
        "fileGaussianInput"    :   "job.inp"      ,
        "fileGaussianOutput"   :   "job.log"      ,
        "pathMopac"            :   os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
        "pathGaussian"         :   os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03")    ,
    }

    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems (): setattr (self, key, value)


    def Run (self):
        # . Read atoms.inp from Molaris
        molaris = MolarisAtomsFile (filename=self.fileAtoms)

        # . Write geometries for the purpose of viewing
        if self.fileTrajectory:
            molaris.WriteQM (filename=self.fileTrajectory, link=True, append=True)

        if self.program == "mopac":
            # . Prepare a MOPAC calculation
            molaris.WriteMopacInput (filename=self.fileMopacInput, method=self.method, charge=self.charge, cosmo=self.cosmo)

            # . Run the calculation
            fileError  = open (self.fileMopacError, "w")
            subprocess.check_call ([self.pathMopac, self.fileMopacInput], stdout=fileError, stderr=fileError)
            fileError.close ()

            # . Convert the output file from MOPAC to forces.out
            mopac = MopacOutputFile (filename=self.fileMopacOutput)
            mopac.WriteMolarisForces (filename=self.fileForces)


        elif self.program == "gaussian":
            # . Prepare a Gaussian calculation
            molaris.WriteGaussianInput (filename=self.fileGaussianInput, methodBasis=self.method, charge=self.charge)

            # . Run the calculation
            fileError  = open (self.fileGaussianError, "w")
            subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError)
            fileError.close ()

            # . Convert the output file from Gaussian to forces.out
            gaussian = GaussianOutputFile (filename=self.fileGaussianOutput)
            gaussian.WriteMolarisForces (filename=self.fileForces)
        else:
            raise exceptions.StandardError ("Program %s not supported." % self.program)
