#-------------------------------------------------------------------------------
# . File      : QMCallerTeraChem.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# . TODO: Electrostatic embedding
#         COSMO model!

import subprocess, os, exceptions, collections, math

from MolarisTools.Parser     import TeraChemOutputFile
from MolarisTools.QMMM       import QMCaller, CS_MULLIKEN, CS_CHELPG, CS_MERZKOLLMAN


Force = collections.namedtuple ("Force", "x  y  z")

_DEFAULT_METHOD = "B3LYP/6-31G*"
_DEFAULT_CONVERGENCE = 5

class QMCallerTeraChem (QMCaller):
    """A class to provide communication between Molaris and TeraChem."""

    defaultAttributes = {
        "restart"                 :   False                 ,
        "method"                  :   _DEFAULT_METHOD       ,
        "SCFConvergence"          :   _DEFAULT_CONVERGENCE  ,
        "fileTeraChemInput"       :   "tc.sp"               ,
        "fileTeraChemOutput"      :   "tc.out"              ,
        "fileTeraChemCheckpoint"  :   "tc.chk"              ,
        "fileCoordinates"         :   "coor.xyz"            ,
        "filePointCharges"        :   "pc.xyz"              ,
        "pathTeraChem"            :   os.path.join (os.environ["HOME"], "TeraChem", "bin", "terachem") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerTeraChem, self).__init__ (**keywordArguments)

        # . Prepare a TeraChem input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write a TeraChem input file."""
        # . Write a coordinates file
        atoms = self.molaris.qatoms + self.molaris.latoms
        natoms = len (atoms)

        output = open (self.fileCoordinates, "w")
        output.write ("%d\nTeraChem job.\n" % natoms)
        for atom in atoms:
            output.write ("%2s    %16.10f    %16.10f    %16.10f\n" % (atom.label, atom.x, atom.y, atom.z))
        output.write ("\n")
        output.close ()

        # . Write an input file
        output = open (self.fileTeraChemInput, "w")
        (method, basis) = self.method.split ("/")
        output.write ("method %s\n" % method)
        output.write ("basis %s\n" % basis)
        output.write ("charge %d\n" % self.charge)
        output.write ("spinmult %d\n" % self.multiplicity)
        output.write ("coordinates %s\n" % self.fileCoordinates)
        output.write ("dftd %s\n" % "no")
        output.write ("run %s\n" % "gradient")

        # . Choose a charge scheme
        schemes = {
            CS_MULLIKEN     :   ""         ,
            CS_MERZKOLLMAN  :   "resp yes" ,
            }
        if (not schemes.has_key (self.chargeScheme)):
            raise exceptions.StandardError ("Charge scheme %s is undefined." % self.chargeScheme)
        selectScheme = schemes[self.chargeScheme]
        if (selectScheme != ""):
            output.write ("%s\n" % selectScheme)

        if (self.qmmm):
            # . Include point charges (electrostatic embedding)
            output.write ("pointcharges %s\n" % self.filePointCharges)

            pcfile = open (self.filePointCharges, "w")
            pointCharges = self.molaris.patoms + self.molaris.watoms
            for atom in pointCharges:
                pcfile.write ("%10.4f  %14.8f  %14.8f  %14.8f\n" % (atom.charge, atom.x, atom.y, atom.z))
            pcfile.close ()
        else:
            # . Do not include point charges
            output.write ("pointcharges %s\n" % "no")

        if (self.fileTeraChemCheckpoint):
            # . Write a checkpoint file
            output.write ("chkfile %s\n" % self.fileTeraChemCheckpoint)

        guess = "generate"
        if (self.restart):
            fileGuess = os.path.join ("scr", "c0")
            if (os.path.exists (fileGuess)):
                guess = fileGuess
            else:
                filea = os.path.join ("scr", "ca")
                fileb = os.path.join ("scr", "cb")
                if (os.path.exists (filea) and os.path.exists (fileb)):
                    guess = "%s %s" % (filea, fileb)
        output.write ("guess %s\n" % guess)

        if (self.SCFConvergence != _DEFAULT_CONVERGENCE):
            output.write ("convthre %e\n" % math.pow (10, -self.SCFConvergence))

        output.write ("end\n\n")
        output.close ()


    def Run (self):
        """Run the calculation."""
        fileOutput = open (self.fileTeraChemOutput, "w")
        subprocess.check_call ([self.pathTeraChem, self.fileTeraChemInput], stdout=fileOutput, stderr=fileOutput)
        fileOutput.close ()
        # . Parse the output file
        terachem = TeraChemOutputFile (filename=self.fileTeraChemOutput)
        self.Efinal = terachem.Efinal

        # . Include forces on QM atoms
        self.forces = terachem.forces

        # . Include forces on point charges
        # if (hasattr (terachem, "pointCharges")):
        #     mmforces = []
        #     for pc in terachem.pointCharges:
        #         force = Force (
        #             x   =   pc.ex   *   pc.charge   ,
        #             y   =   pc.ey   *   pc.charge   ,
        #             z   =   pc.ez   *   pc.charge   , )
        #         mmforces.append (force)
        #     self.mmforces = mmforces

        # . Include charges
        scheme = {
            CS_MULLIKEN     :   terachem.charges    if hasattr (terachem, "charges"   ) else []  ,
            CS_MERZKOLLMAN  :   terachem.espcharges if hasattr (terachem, "espcharges") else []  , }
        self.charges = scheme[self.chargeScheme]
        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
