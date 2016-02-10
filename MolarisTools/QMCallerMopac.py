#-------------------------------------------------------------------------------
# . File      : QMCallerMopac.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from QMCaller            import QMCaller
from Utilities           import WriteData
from MopacOutputFile     import MopacOutputFile

import subprocess, os.path, exceptions


class QMCallerMopac (QMCaller):
    """A class to provide communication between Molaris and Mopac."""

    # . Options specific to Mopac
    # . fileAtoms should be set to "mol.in" if qmmm=True
    defaultAttributes = {
        "fileMopacError"       :   "run.err"      ,
        "fileMopacInput"       :   "run.mop"      ,
        "fileMopacOutput"      :   "run.out"      ,
        "pathMopac"            :   os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerMopac, self).__init__ (**keywordArguments)
        if self.qmmm:
            if self.fileAtoms != "mol.in":
                raise exceptions.StandardError ("With qmmm option enabled, fileAtoms can only be mol.in.")

        # . Prepare a MOPAC input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write a Mopac input file."""
        # . Set up a charge scheme
        schemes = {
            "Mulliken"     :   "MULLIK" ,
            "MerzKollman"  :   "ESP"    ,
            }
        if not schemes.has_key (self.chargeScheme):
            raise exceptions.StandardError ("Charge scheme %s is undefined." % self.chargeScheme)
        chargeScheme = schemes[self.chargeScheme]
        # . Set up a spin state
        multp = {
            1  :   ""        ,
            2  :   "DOUBLET" ,
            3  :   "TRIPLET" ,
            }
        spinState = multp[self.multiplicity]
        # . Set up a solvation model
        solvationModel = ("EPS=%.2f" % self.dielectric) if self.cosmo else ""
        # . Set up a QM/MM model
        qmmmModel = "DEBUG MOL_QMMM" if self.qmmm else ""
        # . Write header
        data   = []
        data.append ("%s  1SCF  CHARGE=%-2d  %s  %s  GRAD  XYZ  %s  %s\n" % (self.method, self.charge, spinState, solvationModel, chargeScheme, qmmmModel))
        data.append ("Comment line\n")
        data.append ("\n")
        # . Write geometry
        atoms = self.molaris.qatoms + self.molaris.latoms
        for atom in atoms:
            data.append ("%2s    %9.4f  1    %9.4f  1    %9.4f  1\n" % (atom.symbol, atom.x, atom.y, atom.z))
        # . Finish up
        WriteData (data, self.fileMopacInput)


    def Run (self):
        """Run the calculation."""
        fileError  = open (self.fileMopacError, "w")
        subprocess.check_call ([self.pathMopac, self.fileMopacInput], stdout=fileError, stderr=fileError)
        fileError.close ()

        # . Parse the output file
        mopac        = MopacOutputFile (filename=self.fileMopacOutput)
        self.Efinal  = mopac.Efinal
        self.forces  = mopac.forces
        if self.chargeScheme == "Mulliken":
            self.charges = mopac.charges
        elif self.chargeScheme == "MerzKollman":
            self.charges = mopac.mkcharges

        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
