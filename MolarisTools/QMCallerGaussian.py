#-------------------------------------------------------------------------------
# . File      : QMCallerGaussian.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from QMCaller            import QMCaller
from Utilities           import WriteData
from GaussianOutputFile  import GaussianOutputFile

import subprocess, os.path, exceptions


class QMCallerGaussian (QMCaller):
    """A class to provide communication between Molaris and Gaussian."""

    # . Options specific to Gaussian
    # . Memory is given in GB
    # . env may define variables such as GAUSS_EXEDIR and GAUSS_SCRDIR
    # . restart means to reuse the wavefunction from the checkpoint file
    defaultAttributes = {
        "env"                     :   None         ,
        "ncpu"                    :   1            ,
        "memory"                  :   1            ,
        "restart"                 :   False        ,
        "fileGaussianError"       :   "job.err"    ,
        "fileGaussianInput"       :   "job.inp"    ,
        "fileGaussianOutput"      :   "job.log"    ,
        "fileGaussianCheckpoint"  :   "job.chk"    ,
        "pathGaussian"            :   os.path.join (os.environ["HOME"], "local", "opt", "g03", "g03") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerGaussian, self).__init__ (**keywordArguments)

        # . Determine if a semiempirical potential is used
        method = self.method[:3]
        if method in ("AM1", "PM3", ) and self.qmmm:
            raise exceptions.StandardError ("Point charges cannot be used with semiempirical methods.")

        # . Reuse the wavefunction if the checkpoint file exists
        if self.fileGaussianCheckpoint:
            self.restart = os.path.exists (self.fileGaussianCheckpoint) and self.restart
        else:
            self.restart = False

        # . Prepare a Gaussian input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write a Gaussian input file."""
        # . Write job control
        data = []
        if self.ncpu > 1:
            data.append ("%%NProcShared=%d\n" % self.ncpu)
        if self.memory:
            data.append ("%%mem=%dgb\n" % self.memory)
        if self.fileGaussianCheckpoint:
            data.append ("%%chk=%s\n" % self.fileGaussianCheckpoint)

        # . Write header
        qmmm    = "CHARGE=ANGSTROMS"          if self.qmmm    else ""
        cosmo   = "SCRF=(Solvent=Water,Read)" if self.cosmo   else ""
        restart = "GUESS=READ"                if self.restart else ""
        data.append ("# %s %s %s %s FORCE NOSYMM\n\n" % (self.method, qmmm, cosmo, restart))
        data.append ("Comment line\n\n")
        data.append ("%d %d\n" % (self.charge, self.multiplicity))

        # . Write geometry
        atoms = self.molaris.qatoms + self.molaris.latoms
        for atom in atoms:
            data.append ("%2s    %9.4f    %9.4f    %9.4f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        data.append ("\n")

        # . If cosmo=True, write epsilon
        if self.cosmo:
            data.append ("eps=%f\n\n" % self.dielectric)

        # . Write point charges
        if self.qmmm:
            pointCharges = self.molaris.patoms + self.molaris.watoms
            for atom in pointCharges:
                data.append ("%9.4f    %9.4f    %9.4f    %9.4f\n" % (atom.x, atom.y, atom.z, atom.charge))
            data.append ("\n")

        # . Finish up
        WriteData (data, self.fileGaussianInput)


    def Run (self):
        # . Run the calculation
        fileError  = open (self.fileGaussianError, "w")
        if self.env:
            subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError, env=self.env)
        else:
            subprocess.check_call ([self.pathGaussian, self.fileGaussianInput], stdout=fileError, stderr=fileError)
        fileError.close ()

        # . Parse the output file
        gaussian     = GaussianOutputFile (filename=self.fileGaussianOutput)
        self.Efinal  = gaussian.Efinal
        self.forces  = gaussian.forces
        self.charges = gaussian.charges

        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
