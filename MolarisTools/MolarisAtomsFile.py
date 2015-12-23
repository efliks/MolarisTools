#-------------------------------------------------------------------------------
# . File      : MolarisAtomsFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units             import *
from   Utilities         import TokenizeLine, WriteData
from   GaussianInputFile import GaussianInputFile
import collections

Atom  = collections.namedtuple ("Atom"  , "symbol x y z charge")
Force = collections.namedtuple ("Force" , "x y z")


class MolarisAtomsFile (object):
    """A class representing atoms for the QC/MM calculation."""

    def __init__ (self, filename="atoms.inp", filterSymbols=["MG", ]):
        """Constructor."""
        self.inputfile     = filename
        self.filterSymbols = filterSymbols
        self._Parse ()


    def _LineToAtom (self, line, includeCharge=False):
        tokens = TokenizeLine (line, converters=[None, float, float, float, float])
        symbol = tokens[0]
        # . Truncate all symbols to one character, except the filtered ones
        atom   = Atom (
            symbol = (symbol[:1] if symbol not in self.filterSymbols else symbol) ,
            charge = (tokens[4] if includeCharge else None) ,
            x      =  tokens[1]  ,
            y      =  tokens[2]  ,
            z      =  tokens[3]  ,)
        return atom


    def _ReadAtoms (self, openfile, natoms, includeCharge=False):
        atoms = []
        for nq in range (natoms):
            line = next (openfile)
            atom = self._LineToAtom (line, includeCharge)
            atoms.append (atom)
        return atoms


    def _Parse (self):
        lines   = open (self.inputfile)
        line    = next (lines)
        # . Get step number and total energy
        step, G = TokenizeLine (line, converters=[int, float])
        try:
            while True:
                line = next (lines)
                # . Read the QM section
                if line.count ("# of qmmm atoms"):
                    nquantum, nlink = TokenizeLine (line, converters=[int, int])
                    # . Read QM atoms proper
                    self.qatoms = self._ReadAtoms (lines, nquantum)
                    # . Read QM link atoms
                    self.latoms = self._ReadAtoms (lines, nlink)
                elif line.count ("# of total frozen protein atoms, # of groups in Region I`"):
                    pass
                elif line.count ("# of frozen water atoms in Region I`"):
                    pass
                # . Read the protein section
                elif line.count ("# of non-frozen protein atoms in Region II"):
                    nprot = TokenizeLine (line, converters=[int])[0]
                    self.patoms = self._ReadAtoms (lines, nprot, includeCharge=True)
                # . Read the free water section
                elif line.count ("# of non-frozen water atoms in the system"):
                    nwater = TokenizeLine (line, converters=[int])[0]
                    self.watoms = self._ReadAtoms (lines, nwater, includeCharge=True)
        except StopIteration:
            pass


    def WriteQM (self, filename="qm.xyz", link=False, append=False):
        """Write QM atoms to an XYZ file."""
        self._WriteSection (filename, (self.qatoms + self.latoms) if link else self.qatoms, includeCharge=False, append=append)


    def WriteProtein (self, filename="prot.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.patoms, includeCharge=True)


    def WriteWater (self, filename="wat.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.watoms, includeCharge=True)


    def _WriteSection (self, filename, atoms, includeCharge=False, skipSymbol=False, append=False):
        """Write an XYZ file."""
        natoms = len (atoms)
        data   = ["%d\n%s\n" % (natoms, filename)]
        for atom in atoms:
            data.append ("%2s    %8.3f    %8.3f    %8.3f    %8s\n" % ("" if skipSymbol else atom.symbol, atom.x, atom.y, atom.z, ("%8.4f" % atom.charge) if includeCharge else ""))
        WriteData (data, filename, append=append)


    def WriteMopacInput (self, filename="run.mop", method="PM3", charge=0, multiplicity=1, eps=78.4, cosmo=False, qmmm=False):
        """Write an input file for MOPAC."""
        multp = {1  :   ""        ,
                 2  :   "DOUBLET" ,
                 3  :   "TRIPLET" }
        data   = []
        data.append ("%s  1SCF  CHARGE=%-2d  %s  %s  GRAD  XYZ  MULLIK  %s\n" % (method, charge, multp[multiplicity], (("EPS=%.2f" % eps) if cosmo else ""), "DEBUG MOL_QMMM" if qmmm else ""))
        data.append ("Comment line\n")
        data.append ("\n")
        for atom in (self.qatoms + self.latoms):
            data.append ("%2s    %9.4f  1    %9.4f  1    %9.4f  1\n" % (atom.symbol, atom.x, atom.y, atom.z))
        WriteData (data, filename)


    def WriteGaussianInput  (self, filename="run.inp", methodBasis="PM3", charge=0, multiplicity=1, qmmm=False):
        """Write an input file for Gaussian."""
        gaussian = GaussianInputFile (
            qmmm           = qmmm          ,
            method         = methodBasis   ,
            charge         = charge        ,
            multiplicity   = multiplicity  ,
            fileInput      = filename      ,)
        gaussian.Write (self.qatoms + self.latoms, pointCharges=(self.patoms + self.watoms))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
