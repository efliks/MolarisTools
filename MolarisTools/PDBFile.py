#-------------------------------------------------------------------------------
# . File      : PDBFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections

PDBChain   = collections.namedtuple ("PDBChain"   , "label  residues")
PDBResidue = collections.namedtuple ("PDBResidue" , "label  serial  chain  atoms")
PDBAtom    = collections.namedtuple ("PDBAtom"    , "label  serial  x  y  z")


class PDBFile (object):
    """A class to read a PDB file."""

    def __init__ (self, filename):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines    = open (self.inputfile)
        # . Initialize
        atoms          = []
        residues       = []
        prevResLabel   = ""
        prevChainLabel = ""
        prevResSerial  = ""
        try:
            while True:
                line = next (lines)
                # . Taken from: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
                #    1 -  6        Record name     "ATOM  "
                #    7 - 11        Integer         Atom serial number.
                #   13 - 16        Atom            Atom name.
                #   17             Character       Alternate location indicator.
                #   18 - 20        Residue name    Residue name.
                #   22             Character       Chain identifier.
                #   23 - 26        Integer         Residue sequence number.
                #   27             AChar           Code for insertion of residues.
                #   31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
                #   39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
                #   47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
                #   55 - 60        Real(6.2)       Occupancy.
                #   61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
                #   73 - 76        LString(4)      Segment identifier, left-justified.
                #   77 - 78        LString(2)      Element symbol, right-justified.
                #   79 - 80        LString(2)      Charge on the atom.
                if line.startswith ("ATOM") or line.startswith ("HETATM"):
                    resLabel   =      line[17:20].strip ()
                    chainLabel =      line[21   ]
                    resSerial  = int (line[22:26])
                    if prevResSerial != resSerial:
                        if prevResSerial != "":
                            residue = PDBResidue (
                                label   =   prevResLabel    ,
                                serial  =   prevResSerial   ,
                                chain   =   prevChainLabel  ,
                                atoms   =   atoms           ,
                                )
                            residues.append (residue)
                            atoms = []
                    atom = PDBAtom (
                        label   =          line[12:16].strip () ,
                        serial  =     int (line[6 :11]) ,
                        x       =   float (line[30:38]) ,
                        y       =   float (line[38:46]) ,
                        z       =   float (line[46:54]) ,
                        )
                    atoms.append (atom)
                    prevResLabel   = resLabel
                    prevChainLabel = chainLabel
                    prevResSerial  = resSerial
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . Finish up
        if atoms:
            residue = PDBResidue (
                label   =   prevResLabel    ,
                serial  =   prevResSerial   ,
                chain   =   prevChainLabel  ,
                atoms   =   atoms           ,
                )
            residues.append (residue)
        self.residues = residues


    def WriteResidue (self, resSerial, filename="residue.pdb", segLabel="PRTA"):
        """Write a residue to a separate PDB file."""
        output = []
        for residue in self.residues:
            if residue.serial == resSerial:
                for atom in residue.atoms:
                    # . -3s replace by -4s and remove one heading space
                    output.append ("%-6s%5d  %-3s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%22s\n" % (
                        "ATOM", atom.serial, 
                        atom.label, "", residue.label, 
                        residue.chain, residue.serial, "", 
                        atom.x, atom.y, atom.z, segLabel))
                output.append ("TER\n")
                output.append ("END\n")
                break
        if output:
            fo = open (filename, "w")
            for line in output:
                fo.write (line)
            fo.close ()


    @property
    def nresidues (self):
        if hasattr (self, "residues"):
            return len (self.residues)
        return 0

    @property
    def natoms (self):
        if hasattr (self, "residues"):
            natoms = 0
            for residue in self.residues:
                natoms += len (residue.atoms)
            return natoms
        return 0


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
