#-------------------------------------------------------------------------------
# . File      : QMCallerQChem.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from QMCaller         import QMCaller, CS_MULLIKEN, CS_CHELPG, CS_MERZKOLLMAN
from QChemOutputFile  import QChemOutputFile, EfieldFile

import subprocess, os.path, exceptions, collections


Force = collections.namedtuple ("Force", "x  y  z")

_DEFAULT_SAV_FOLDER = "sav"


class QMCallerQChem (QMCaller):
    """A class to provide communication between Molaris and Q-Chem."""

    # . Options specific to Q-Chem
    defaultAttributes = {
        "job"           :   "job"     ,
        "ncpu"          :   1         ,
        "memory"        :   1         ,
        "restart"       :   False     ,
        "exchange"      :   ""        ,
        "correlation"   :   ""        ,
        "basis"         :   ""        ,
        "method"        :   "B3LYP/6-31G*" ,
        "scratch"       :   os.path.join (os.environ["QCSCRATCH"]) if os.environ.has_key ("QCSCRATCH") else "." ,
        "pathQChem"     :   os.path.join (os.environ["QC"])        if os.environ.has_key ("QC")        else ""  ,
            }
    defaultAttributes.update (QMCaller.defaultAttributes)


    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (QMCallerQChem, self).__init__ (**keywordArguments)
        # . Prepare a Q-Chem input file
        self._WriteInput ()


    def _WriteInput (self):
        """Write QChem input files."""
        # . Check for scratch space
        if not os.path.exists (self.scratch):
            os.makedirs (self.scratch)
        # . Header
        lines = ["$comment", ]
        lines.append ("Q-Chem Job.")
        lines.append ("$end")
        lines.append ("")
        # . Control section
        lines.append ("$rem")
        lines.append ("JOBTYPE FORCE")
        if self.qmmm:
            lines.append ("QM_MM TRUE")
            lines.append ("QMMM_PRINT TRUE")
        if self.exchange != "":
            lines.append ("EXCHANGE %s" % self.exchange)
        if self.correlation != "":
            lines.append ("CORRELATION %s" % self.correlation)
        if self.basis != "":
            lines.append ("BASIS %s" % self.basis)
        if self.method != "":
            (method, basis) = self.method.split ("/")
            lines.append ("METHOD %s" % method)
            lines.append ("BASIS %s" % basis)
        lines.append ("MEM_TOTAL %d" % (self.memory * 1000))
        lines.append ("SCF_CONVERGENCE 5")
        lines.append ("THRESH 12")
        lines.append ("SYMMETRY OFF")
        lines.append ("SYM_IGNORE TRUE")
        if self.restart:
            if os.path.exists (os.path.join (self.scratch, _DEFAULT_SAV_FOLDER)):
                lines.append ("SCF_GUESS READ")
        lines.append ("$end")
        # . Molecule geometry
        lines.append ("")
        lines.append ("$molecule")
        lines.append ("%d %d" % (self.charge, self.multiplicity))
        atoms = self.molaris.qatoms + self.molaris.latoms
        for atom in atoms:
            lines.append ("%2s    %16.10f    %16.10f    %16.10f" % (atom.label, atom.x, atom.y, atom.z))
        lines.append ("$end")
        lines.append ("")
        # . Point charges
        if self.qmmm:
            lines.append ("$external_charges")
            pointCharges = self.molaris.patoms + self.molaris.watoms
            for atom in pointCharges:
                lines.append ("%16.10f    %16.10f    %16.10f    %12.4f" % (atom.x, atom.y, atom.z, atom.charge))
            lines.append ("$end")
        lines.append ("")
        # . Write everything to a file
        fo = open (os.path.join (self.scratch, (self.job + ".inp")), "w")
        for line in lines:
            fo.write (line + "\n")
        fo.close ()


    def Run (self):
        """Run the calculation."""
        qchemInput  =  os.path.join (self.scratch  ,  self.job + ".inp")
        qchemOutput =  os.path.join (self.scratch  ,  self.job + ".out")
        qchemError  =  os.path.join (self.scratch  ,  self.job + ".err")
        qchemEField =  os.path.join (self.scratch  ,  "efield.dat"     )
        # . Call Q-Chem
        fileError   = open (qchemError, "w")
        if self.ncpu < 2:
            command = [os.path.join (self.pathQChem, "bin", "qchem"), "-save", qchemInput, qchemOutput, _DEFAULT_SAV_FOLDER]
        else:
            command = [os.path.join (self.pathQChem, "bin", "qchem"), "-save", "-nt", "%d" % self.ncpu, qchemInput, qchemOutput, _DEFAULT_SAV_FOLDER]
        subprocess.check_call (command, stdout=fileError, stderr=fileError)
        fileError.close ()
        # . Parse output files
        qchem   = QChemOutputFile (qchemOutput)
        efield  = EfieldFile (qchemEField)
        if self.qmmm:
            # . Calculate electrostatic forces acting on MM atoms
            mmforces     = []
            pointCharges = self.molaris.patoms + self.molaris.watoms
            nvectors     = len (pointCharges)
            for point, (ex, ey, ez) in zip (pointCharges, efield.field[:nvectors]):
                force = Force (
                    x   =   -ex  *  point.charge   ,
                    y   =   -ey  *  point.charge   ,
                    z   =   -ez  *  point.charge   , )
                mmforces.append (force)
            self.mmforces = mmforces
        # . Include forces on QM atoms
        forces = []
        for (fx, fy, fz) in efield.field[nvectors:]:
            force = Force (
                x  =  -fx  ,
                y  =  -fy  ,
                z  =  -fz  , )
            forces.append (force)
        self.forces = forces
        # . If there are point charges, remove their self interaction energy from the final QM energy
        self.Efinal  = (qchem.Efinal - qchem.Echrg) if self.qmmm else qchem.Efinal
        # . Include charges
        if   self.chargeScheme == CS_MULLIKEN:
            self.charges = qchem.charges
        elif self.chargeScheme == CS_MERZKOLLMAN:
            raise exceptions.StandardError ("Merz-Kollman charges are not (yet) implemented in QMCallerQChem.")
        elif self.chargeScheme == CS_CHELPG:
            raise exceptions.StandardError ("CHELPG charges are not (yet) implemented in QMCallerQChem.")
        # . Finish up
        self._Finalize ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
