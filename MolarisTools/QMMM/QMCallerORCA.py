#-------------------------------------------------------------------------------
# . File      : QMCallerORCA.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import subprocess, os.path, exceptions

from MolarisTools.Utilities  import WriteData
from MolarisTools.Parser     import ORCAOutputFile, PCgradFile, EngradFile
from MolarisTools.QMMM       import QMCaller, CS_MULLIKEN, CS_CHELPG, CS_MERZKOLLMAN


class QMCallerORCA (QMCaller):
    """A class to provide communication between Molaris and ORCA."""

    # . Options specific to ORCA.
    # . Note that ORCA will by default reuse the previous orbitals as a guess for SCF, hence no restart option.
    defaultAttributes = {
        "job"        :   "job"          ,
        "scratch"    :   "orca"         ,
        "ncpu"       :   1              ,
        "memory"     :   1              ,
        "method"     :   "B3LYP/6-31G*" ,
        "debug"      :   False          ,
        "pathORCA"   :   os.path.join (os.environ["HOME"], "local", "opt", "orca_3_0_0_linux_x86-64", "orca") ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerORCA, self).__init__ (**keywordArguments)
        # . Prepare a ORCA input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write ORCA input files."""
        # . Check for the scratch directory
        if not os.path.exists (self.scratch):
            os.makedirs (self.scratch)
        # . Header
        lines  = ["# . ORCA job", ]
        # . Include solvent or protein
        if   self.qmmm:
            pcFile = os.path.abspath (os.path.join (self.scratch, self.job + ".pc"))
            lines.append ("%%pointcharges \"%s\"" % pcFile)
        elif self.cosmo:
            raise exceptions.StandardError ("COSMO model is not (yet) implemented in QMCallerORCA.")
        # . Number of processors
        if self.ncpu < 2:
            cpus = ""
        else:
            cpus = " PAL%d" % self.ncpu
        # . Level of theory
        method, basis = self.method.split ("/")
        lines.append ("! ENGRAD %s %s SCFCONV10%s" % (method, basis, cpus))
        # . Electronic state
        lines.append ("* xyz %d %d" % (self.charge, self.multiplicity))
        # . Geometry
        atoms = self.molaris.qatoms + self.molaris.latoms
        for atom in atoms:
            lines.append ("%2s    %16.10f    %16.10f    %16.10f" % (atom.label, atom.x, atom.y, atom.z))
        # . End of file
        lines.append ("*")
        # . Write everything to a file
        fo = open (os.path.join (self.scratch, (self.job + ".inp")), "w")
        for line in lines:
            fo.write (line + "\n")
        fo.close ()
        # . Now prepare PC data
        if self.qmmm:
            pointCharges = self.molaris.patoms + self.molaris.watoms
            ncharges = len (pointCharges)
            lines    = ["  %d" % ncharges, ]
            for atom in pointCharges:
                lines.append ("%12.4f    %16.10f    %16.10f    %16.10f" % (atom.charge, atom.x, atom.y, atom.z))
            # . Write point charges to a file
            fo = open (os.path.join (self.scratch, (self.job + ".pc")), "w")
            for line in lines:
                fo.write (line + "\n")
            fo.close ()


    def Run (self):
        # . Run the calculation
        orcaInput   =   os.path.join (self.scratch, self.job + ".inp")
        orcaOutput  =   os.path.join (self.scratch, self.job + ".log")
        # . In the debug mode, reuse the already existing log file
        calculate   = True
        if self.debug:
            if os.path.exists (orcaOutput):
                calculate = False
        if calculate:
            fileOutput  = open (orcaOutput, "w")
            subprocess.check_call ([self.pathORCA, orcaInput], stdout=fileOutput, stderr=fileOutput)
            fileOutput.close ()
        # . Parse the output file
        orca        = ORCAOutputFile (orcaOutput, reverse=True)
        # . In ORCA, the final QM energy does not seem to include the self interaction energy of point charges
        self.Efinal = orca.Efinal
        # . Include forces on QM atoms
        engrad      = EngradFile (os.path.join (self.scratch, self.job + ".engrad"), reverse=True)
        self.forces = engrad.forces
        # . Include forces on point charges
        if self.qmmm:
            pcgrad        = PCgradFile (os.path.join (self.scratch, self.job + ".pcgrad"), reverse=True)
            self.mmforces = pcgrad.forces
        # . Include charges
        if   self.chargeScheme == CS_MULLIKEN:
            charges = []
            for atom in orca.atoms:
                charges.append (atom.charge)
            self.charges = charges
        elif self.chargeScheme == CS_MERZKOLLMAN:
            raise exceptions.StandardError ("Merz-Kollman charges are not (yet) implemented in QMCallerORCA.")
        elif self.chargeScheme == CS_CHELPG:
            raise exceptions.StandardError ("CHELPG charges are not (yet) implemented in QMCallerORCA.")
        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
