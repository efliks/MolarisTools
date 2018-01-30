#-------------------------------------------------------------------------------
# . File      : TeraChemOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections, exceptions, os

from  MolarisTools.Units     import HARTREE_TO_KCAL_MOL, HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, GRADIENT_TO_FORCE
from  MolarisTools.Utilities import TokenizeLine

Atom = collections.namedtuple ("Atom", "symbol  x  y  z  charge")
Force = collections.namedtuple ("Force", "x  y  z")


class TeraChemOutputFile (object):
    """A class to read a TeraChem output file."""

    def __init__ (self, filename="tc.out", deep=True):
        """Constructor."""
        self.inputfile = filename
        self.deep = deep
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                if line.startswith ("FINAL ENERGY:"):
                    tokens = TokenizeLine (line, converters=[None, None, float, None])
                    self.Efinal = tokens[2] * HARTREE_TO_KCAL_MOL

                elif line.startswith ("Scratch directory:"):
                    tokens = TokenizeLine (line, converters=[None, None, None])
                    self.scratchFolder = tokens[2]

                elif line.startswith ("Total atoms:"):
                    tokens = TokenizeLine (line, converters=[None, None, int])
                    natoms = tokens[2]
                    # if (natoms != len (tempatoms)):
                    #     pass

                # ****** QM coordinates ******
                # C        -0.0178447840         0.0103903440        -0.0015978260  
                # H        -0.0463346130         0.0459578690         1.0997854660  
                # H         1.0186858510         0.0532897140        -0.3692509610  
                # H        -0.5543873770        -0.8695391090        -0.3892902580  
                # Cl        -0.8673247020         1.4796754130        -0.6190902640  
                # Cl         1.5376616510        -2.4337214350         0.8399845050  
                # 
                elif line.startswith ("****** QM coordinates ******"):
                    tempatoms = []
                    while True:
                        line = next (lines)
                        if (line.strip () == ""):
                            break
                        tokens = TokenizeLine (line, converters=[None, float, float, float])
                        atom = Atom (symbol=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], charge=0.0)
                        tempatoms.append (atom)

                # ESP unrestraint charges:
                # Atom          X          Y          Z     Charge   Exposure
                # -----------------------------------------------------------
                #    C  -0.033722   0.019635  -0.003019  -0.161628     0.0529
                #    H  -0.087560   0.086848   2.078293   0.135770     0.5000
                #    H   1.925037   0.100703  -0.697783   0.150080     0.4130
                #    H  -1.047640  -1.643191  -0.735652   0.144291     0.4348
                #   Cl  -1.639006   2.796181  -1.169911  -0.328851     0.8444
                #   Cl   2.905759  -4.599067   1.587341  -0.939662     0.9270
                # -----------------------------------------------------------
                elif line.startswith ("ESP unrestraint charges:"):
                    next (lines)
                    next (lines)
                    self.espcharges = []
                    while True:
                        line = next (lines)
                        if line.startswith ("----"):
                            break
                        tokens = TokenizeLine (line, converters=[None, float, float, float, float, float])
                        (symbol, charge) = (tokens[0], tokens[4])
                        self.espcharges.append (charge)

                # Gradient units are Hartree/Bohr
                # ---------------------------------------------------
                #         dE/dX            dE/dY            dE/dZ
                #    -0.0094514443     0.0141415195    -0.0068751376
                #    -0.0018028819     0.0035487049     0.0113054990
                #     0.0100877721     0.0051162387    -0.0050948764
                #    -0.0088062409    -0.0065296736    -0.0054152317
                #     0.0115697382    -0.0189485242     0.0072745597
                #    -0.0015969432     0.0026717342    -0.0011948132
                # ---------------------------------------------------
                elif line.startswith ("Gradient units are Hartree/Bohr"):
                    next (lines)
                    next (lines)
                    self.forces = []
                    while True:
                        line = next (lines)
                        if (line.startswith ("----") or line.strip () == ""):
                            break
                        tokens = TokenizeLine (line, converters=[float, float, float])
                        (fx, fy, fz) = tokens
                        force = Force (x=(fx * GRADIENT_TO_FORCE * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM), y=(fy * GRADIENT_TO_FORCE * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM), z=(fz * GRADIENT_TO_FORCE * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM))
                        self.forces.append (force)

        except exceptions.StopIteration:
            pass
        lines.close ()

        if (self.deep):
            # . Deep parsing (aka parse files in the scratch folder)
            if hasattr (self, "scratchFolder"):

                # . Collect Mulliken charges
                lines = open (os.path.join (self.scratchFolder, "charge_mull.xls"))
                self.charges = []
                for i in range (natoms):
                    line = next (lines)
                    tokens = TokenizeLine (line, converters=[int, None, float])
                    self.charges.append (tokens[2])
                lines.close ()

                # . Collect coordinates
                fileGeometry = os.path.join (self.scratchFolder, "xyz.xyz")
                if (os.path.exists (fileGeometry)):
                    lines = open (fileGeometry)
                    next (lines)
                    next (lines)
                    self.atoms = []
                    for i in range (natoms):
                        line = next (lines)
                        tokens = TokenizeLine (line, converters=[None, float, float, float])
                        atom = Atom (symbol=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], charge=self.charges[i])
                        self.atoms.append (atom)
                    lines.close ()

        # . Finalize after reading all files
        # if (not hasattr (self, "atoms")):
        #     self.atoms = tempatoms

    @property
    def natoms (self):
        if (hasattr (self, "atoms")):
            return len (self.atoms)
        return 0

    @property
    def ncharges (self):
        if (hasattr (self, "pointCharges")):
            return len (self.pointCharges)
        return 0

    def WriteMolarisForces (self, filename="forces.out", Eref=0., useESPCharges=False):
        """Write a file in the Molaris-suitable format."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
