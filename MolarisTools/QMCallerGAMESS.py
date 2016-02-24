#-------------------------------------------------------------------------------
# . File      : QMCallerGAMESS.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Units               import symbolToAtomicNumber
from QMCaller            import QMCaller
from Utilities           import WriteData
from GAMESSOutputFile    import GAMESSOutputFile, GAMESSDatFile

import subprocess, os.path, exceptions


class QMCallerGAMESS (QMCaller):
    """A class to provide communication between Molaris and GAMESS-US."""

    # . Options specific to GAMESS
    # . By default, basis set is 6-31G*
    defaultAttributes = {
        "ncpu"                  :   1            ,
        "version"               :   "01"         ,
        "memory"                :   10           ,
        "restart"               :   False        ,
        "gbasis"                :   "n31"        ,
        "ngauss"                :   6            ,
        "ndfunc"                :   1            ,
        "fileGAMESSError"       :   "job.err"    ,
        "fileGAMESSInput"       :   "job.inp"    ,
        "fileGAMESSOutput"      :   "job.log"    ,
        "fileGAMESSCheckpoint"  :   "job.dat"    ,
        "pathGAMESS"            :   os.path.join (os.environ["HOME"], "local", "opt", "gamess", "rungms") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerGAMESS, self).__init__ (**keywordArguments)
        # . Prepare a GAMESS input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write a GAMESS input file."""
        data = []
        # . System section (memory is in MW)
        data.append (" $system mwords=%d $end\n" % self.memory)

        # . Control section, charge and multiplicity
        data.append (" $contrl scftyp=rhf runtyp=gradient dfttyp=%s\n" % self.method)
        data.append ("      maxit=100 mult=%d icharg=%d $end\n" % (self.multiplicity, self.charge))

        # . Basis set
        data.append (" $basis  gbasis=%s ngauss=%d ndfunc=%d $end\n" % (self.gbasis, self.ngauss, self.ndfunc))

        # . Initial guess
        # . Reuse the wavefunction if the checkpoint file exists
        guess = " $guess  guess=huckel $end\n"
        if self.restart:
            if os.path.exists (self.fileGAMESSCheckpoint):
                guess = " $guess  guess=moread $end\n"
            else:
                self.restart = False
        data.append (guess)

        # . SCF options
        data.append (" $scf    dirscf=.true. $end\n")

        # . Geometry
        data.append (" $data\n")
        data.append ("Comment line\n")
        data.append ("C1\n")
        atoms = self.molaris.qatoms + self.molaris.latoms
        for atom in atoms:
            data.append ("%2s    %4.1f    %9.4f    %9.4f    %9.4f\n" % (atom.label, symbolToAtomicNumber[atom.label], atom.x, atom.y, atom.z))
        data.append ("$end\n")

        # . Initial orbitals
        if self.restart:
            data.append (" $vec\n")
            checkpoint = GAMESSDatFile (self.fileGAMESSCheckpoint)
            data.extend (checkpoint.vec)
            data.append ("$end\n")

        # . Finish up
        WriteData (data, self.fileGAMESSInput)

        # . TODO: Cosmo and QM/MM (point charges)


    def Run (self):
        # . Run the calculation
        # . Open files
        fileError       = open (self.fileGAMESSError  , "w")
        fileOutput      = open (self.fileGAMESSOutput , "w")
        stem, extension = os.path.splitext (self.fileGAMESSInput)
        filename        = os.path.basename (stem)
        # . Example: ~/local/opt/gamess/rungms scf-ginkgo 01 8 > scf-ginkgo.out &
        subprocess.check_call ([self.pathGAMESS, filename, self.version, "%d" % self.ncpu], stdout=fileOutput, stderr=fileError)
        # . Close files
        fileOutput.close ()
        fileError.close  ()

        # . Parse the output file
        gamess       = GAMESSOutputFile (filename=self.fileGAMESSOutput)
        self.Efinal  = gamess.Efinal
        self.forces  = gamess.forces
        # . Assign charges
        if self.chargeScheme == "Mulliken":
            self.charges = gamess.charges

        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
