#-------------------------------------------------------------------------------
# . File      : MolarisAtomsFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections, math, os

from MolarisTools.Utilities  import TokenizeLine, WriteData


Atom  = collections.namedtuple ("Atom"  , "label  charge  x  y  z")
Force = collections.namedtuple ("Force" , "x  y  z")



class MolarisAtomsFile (object):
    """A class representing atoms for the QC/MM calculation."""


    def CalculatePotentials (self):
        """Calculate the electrostatic potential at centers of the QM atoms."""
        # . Solute QM atoms, including link atoms
        qmatoms = self.qatoms + self.latoms

        # . Protein and waters MM atoms
        mmatoms = self.patoms + self.watoms

        # . Calculate the electrostatic potential for each QM center
        for qm in qmatoms:
            elpot = 0.
            for mm in mmatoms:
                r = math.sqrt ((qm.x - mm.x) ** 2 + (qm.y - mm.y) ** 2 + (qm.z - mm.z) ** 2)
                elpot += (mm.charge / r)
            print ("%2s    %16.10f    %16.10f    %16.10f    %16.10f" % (qm.label, qm.x, qm.y, qm.z, 332. * elpot))


    def __init__ (self, filename="mol.in", replaceSymbols=None):
        """Constructor."""
        self.inputfile      = filename
        self.replaceSymbols = replaceSymbols
        self._Parse ()


    def _LineToAtom (self, line, includeCharge=False):
        # P     4.068890911    1.177766760   16.194718196  157.867317447
        # O     5.406419490    1.617535751   16.790930680  182.964771865
        # ...
        # 
        # tokens = TokenizeLine (line, converters=[None, float, float, float, float])
        # . Assume a column-based format, instead of open one
        tokens = (line[:2], float (line[2:18]), float (line[18:32]), float (line[32:47]), float (line[47:62]))

        symbol = tokens[0]
        # . Replace selected atomic symbols (a workaround for a persisting bug in Molaris)
        if self.replaceSymbols:
            for (symbolOld, symbolNew) in self.replaceSymbols:
                if symbolOld == symbol:
                    symbol = symbolNew
                    break
        # . Create an atom
        atom   = Atom (
            label  =  symbol,
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
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Read the header
                if   line.count ("MD step"):
                    line    = line[:line.find ("MD step")]
                    tokens  = TokenizeLine (line)
                    if len (tokens) != 9:
                        # . Old format
                        tokens       = TokenizeLine (line, converters=[int, float, float, float, float, float, float, int])
                        step, energy = tokens[:2]
                    else:
                        # . New format
                        tokens       = TokenizeLine (line, converters=[int, float, float, float, float, float, float, int, int])
                        (step, energy), stateID = tokens[:2], tokens[-1]
                        self.stateID = stateID
                    self.mdstep  = step
                    self.Etot    = energy

                # . Read the QM section
                elif line.count ("# of qmmm atoms"):
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
        # . Close the file
        lines.close ()


    def WriteQM (self, filename="qm.xyz", link=False, append=False, caption=""):
        """Write QM atoms to an XYZ file."""
        self._WriteSection (filename, (self.qatoms + self.latoms) if link else self.qatoms, includeCharge=False, append=append, caption=caption)


    def WriteProtein (self, filename="prot.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.patoms, includeCharge=True)


    def WriteWater (self, filename="wat.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.watoms, includeCharge=True)


    def _WriteSection (self, filename, atoms, includeCharge=False, skipSymbol=False, append=False, caption=""):
        """Write an XYZ file."""
        natoms = len (atoms)
        data   = ["%d\n%s\n" % (natoms, caption if caption else filename)]
        for atom in atoms:
            data.append ("%2s    %8.3f    %8.3f    %8.3f    %8s\n" % ("" if skipSymbol else atom.label, atom.x, atom.y, atom.z, ("%8.4f" % atom.charge) if includeCharge else ""))
        WriteData (data, filename, append=append)


    def WriteMopacInput (self, filename="run.mop", method="PM3", charge=0, multiplicity=1, eps=78.4, cosmo=False, qmmm=False):
        """Write an input file for MOPAC."""
        pass


    def WriteGaussianInput  (self, filename="run.inp", methodBasis="PM3", charge=0, multiplicity=1, qmmm=False):
        """Write an input file for Gaussian."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"):
    # . Read in a test file
    if (os.path.exists ("mol.in")):
        molin = MolarisAtomsFile ()
        molin.WriteQM ()
