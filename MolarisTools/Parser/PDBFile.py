#-------------------------------------------------------------------------------
# . File      : PDBFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import  collections, exceptions

PDBChain   = collections.namedtuple ("PDBChain"   , "label  residues")
PDBAtom    = collections.namedtuple ("PDBAtom"    , "label  serial  x  y  z")

_MODULE_LABEL      = "PDBFile"
_DEFAULT_LOG_LEVEL = 1
_PDB_FORMAT_ATOM   = "%-6s%5d %-4s %3s %1s%4s%1s   %8.3f%8.3f%8.3f%22s\n"


class PDBResidue (object):
    """A class to represent a PDB residue."""

    # label  serial  chain  atoms  parent
    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in keywordArguments.iteritems ():
            setattr (self, key, value)


    def GetBonds (self):
        """Return a list of bonds for a residue."""
        pairs = []
        pdb   = self.parent
        if pdb.bonds != []:
            for (seriala, serialb) in pdb.bonds:
                found = False
                for i, atoma in enumerate (self.atoms):
                    if atoma.serial == seriala:
                        found = True
                        break
                if not found:
                    continue
                paira = (i, atoma.serial, atoma.label)
                found = False
                for j, atomb in enumerate (self.atoms):
                    if atomb.serial == serialb:
                        found = True
                        break
                if not found:
                    continue
                pairb = (j, atomb.serial, atomb.label)
                pair  = (paira, pairb)
                pairs.append (pair)
        return pairs


    def ReplaceAtom (self, label, newLabel):
        """Replace an atom."""
        found = False
        for (i, atom) in enumerate (self.atoms):
            if atom.label == label:
                new = PDBAtom (
                    label   =   newLabel    ,
                    serial  =   atom.serial ,
                    x       =   atom.x      ,
                    y       =   atom.y      ,
                    z       =   atom.z      , )
                self.atoms[i] = new
                found = True
        if not found:
            raise exceptions.StandardError ("Atom %s not found." % label)


    def KillAtom (self, label):
        """Replace an atom."""
        new   = []
        found = False
        for atom in self.atoms:
            if atom.label == label:
                found = True
                continue
            new.append (atom)
        if not found:
            raise exceptions.StandardError ("Atom %s not found." % label)
        self.atoms = new
        # . Remove bonds that involve removed atom
        remove = []
        for (paira, pairb) in self.bonds:
            if (paira == label) or (pairb == label):
                pair = (paira, pairb)
                remove.append (pair)
        pdb    = self.parent
        update = []
        for (paira, pairb) in pdb.bonds:
            if ((paira, pairb) in remove) or ((pairb, paira) in remove):
                continue
            pair = (paira, pairb)
            update.append (pair)
        pdb.bonds = update


    @property
    def bonds (self):
        return self.GetBonds ()

    @property
    def nbonds (self):
        bonds = self.GetBonds ()
        return len (bonds)


#===============================================================================
class PDBFile (object):
    """A class to read a PDB file."""


    def CenterMolecule (self, residueSerial, atomSerial):
        """Move a selected atom to the origin of the coordinate system."""
        found = False
        for residue in self.residues:
            if (residue.serial == residueSerial):
                for atom in residue.atoms:
                    if (atom.serial == atomSerial):
                        (ax, ay, az) = (atom.x, atom.y, atom.z)
                        found = True
                        break
                break
        if not found:
            raise exceptions.StandardError ("Residue or atom not found.")
        for residue in self.residues:
            newAtoms = []
            for atom in residue.atoms:
                newAtom = PDBAtom (
                    label   =  atom.label     ,
                    serial  =  atom.serial    ,
                    x       =  (atom.x - ax)  ,
                    y       =  (atom.y - ay)  ,
                    z       =  (atom.z - az)  , )
                newAtoms.append (newAtom)
            residue.atoms = newAtoms


    def __init__ (self, filename, logLevel=_DEFAULT_LOG_LEVEL):
        """Constructor."""
        self.inputfile = filename
        self.logLevel  = logLevel
        self._Parse ()


    def _Parse (self):
        lines    = open (self.inputfile)
        if self.logLevel > 1:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.inputfile))
        # . Initialize
        bonds          = []
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
                                parent  =   self            ,
                                )
                            residues.append (residue)
                            atoms = []
                            if self.logLevel > 1:
                                print ("# . %s> Added residue %s %s %d" % (_MODULE_LABEL, prevChainLabel, prevResLabel, prevResSerial))
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


                #  COLUMNS         DATA TYPE        FIELD           DEFINITION
                #  ---------------------------------------------------------------------------------
                #   1 -  6         Record name      "CONECT"
                #  
                #   7 - 11         Integer          serial          Atom serial number
                #  
                #  12 - 16         Integer          serial          Serial number of bonded atom
                #  
                #  17 - 21         Integer          serial          Serial number of bonded atom
                #  
                #  22 - 26         Integer          serial          Serial number of bonded atom
                #  
                #  27 - 31         Integer          serial          Serial number of bonded atom
                #  
                #  32 - 36         Integer          serial          Serial number of hydrogen bonded
                #                                                   atom
                #  
                #  37 - 41         Integer          serial          Serial number of hydrogen bonded
                #                                                   atom
                #  
                #  42 - 46         Integer          serial          Serial number of salt bridged
                #                                                   atom
                #  
                #  47 - 51         Integer          serial          Serial number of hydrogen bonded
                #                                                   atom
                #  
                #  52 - 56         Integer          serial          Serial number of hydrogen bonded
                #                                                   atom
                #  
                #  57 - 61         Integer          serial          Serial number of salt bridged
                #                                                   atom
                elif line.startswith ("CONECT"):
                    token   =  line[ 6:11]
                    serial  =  int (token)
                    tokens  =  line[11:16], line[16:21], line[21:26], line[26:31]
                    serialsBonded = []
                    for token in tokens:
                        token = token.strip ()
                        if token.isdigit ():
                            serialsBonded.append (int (token))
                    if len (serialsBonded) > 0:
                        for serialBonded in serialsBonded:
                            pair    = (serial       , serialBonded)
                            pairRev = (serialBonded , serial      )
                            if pair not in bonds:
                                if pairRev not in bonds:
                                    seriala, serialb = pair
                                    if seriala > serialb:
                                        seriala, serialb = serialb, seriala
                                        pair = (seriala, serialb)
                                    bonds.append (pair)
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
                parent  =   self            ,
                )
            residues.append (residue)
            if self.logLevel > 1:
                print ("# . %s> Added residue %s %s %d" % (_MODULE_LABEL, prevChainLabel, prevResLabel, prevResSerial))
        if self.logLevel > 0:
            nresidues = len (residues)
            print ("# . %s> Found %d residues" % (_MODULE_LABEL, nresidues))
            nbonds    = len (bonds)
            print ("# . %s> Found %d bonds" % (_MODULE_LABEL, nbonds))
        self.residues = residues
        self.bonds    = bonds


    def Write (self, filename="protein.pdb"):
        """Write all residues to one file."""
        output     = []
        atomSerial = 1
        for residue in self.residues:
            (atomSerial, lines) = self.WriteResidue (residue.serial, terminate=False, startSerial=atomSerial)
            output.extend (lines)
        output.append ("TER\n")
        output.append ("END\n")
        if output:
            fo = open (filename, "w")
            for line in output:
                fo.write (line)
            fo.close ()


    def WriteResidue (self, resSerial, filename="residue.pdb", segLabel="PRTA", terminate=True, startSerial=1):
        """Write a residue to a separate PDB file."""
        output = []
        for residue in self.residues:
            if residue.serial == resSerial:
                for iatom, atom in enumerate (residue.atoms):
                    atomLabel  = atom.label if len (atom.label) > 3 else (" %s" % atom.label)
                    atomSerial = iatom + startSerial
                    output.append (_PDB_FORMAT_ATOM % ("ATOM", atomSerial, atomLabel, residue.label, residue.chain, residue.serial, "", atom.x, atom.y, atom.z, segLabel))
                if terminate:
                    output.append ("TER\n")
                    output.append ("END\n")
                break
        if output:
            if terminate:
                fo = open (filename, "w")
                for line in output:
                    fo.write (line)
                fo.close ()
        return (atomSerial + 1, output)


    def CheckForMissingAtoms (self, library, includeHydrogens=False):
        """Check all residues if they are missing any atoms."""
        for residue in self.residues:
            if library.has_key (residue.label):
                # . Collect labels from the library
                component     = library[residue.label]
                libraryLabels = []
                for atom in component.atoms:
                    libraryLabels.append (atom.atomLabel)
                # . Collect labels from PDB atoms
                pdbLabels     = []
                for atom in residue.atoms:
                    pdbLabels.append (atom.label)
                # . Check if all labels are present in the PDB file
                missing   = []
                for libraryLabel in libraryLabels:
                    if not includeHydrogens:
                        if libraryLabel[0] == "H":
                            continue
                    if libraryLabel not in pdbLabels:
                        missing.append (libraryLabel)
                if self.logLevel > 0:
                    if missing:
                        nmissing = len (missing)
                        atoms    = " ".join (missing)
                        print ("Residue %s %s %d has %d missing atom%s: %s" % (residue.chain, residue.label, residue.serial, nmissing, "s" if nmissing > 1 else "", atoms))
                # . Check for redundant labels in the PDB file
                redundant = []
                for pdbLabel in pdbLabels:
                    if not includeHydrogens:
                        if pdbLabel[0] == "H":
                            continue
                    if pdbLabel not in libraryLabels:
                        redundant.append (pdbLabel)
                if self.logLevel > 0:
                    if redundant:
                        nredundant = len (redundant)
                        atoms      = " ".join (redundant)
                        print ("Residue %s %s %d has %d redundant atom%s: %s" % (residue.chain, residue.label, residue.serial, nredundant, "s" if nredundant > 1 else "", atoms))
            else:
                if self.logLevel > 0:
                    print ("Residue %s %s %d not found in the library." % (residue.chain, residue.label, residue.serial))


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

    @property
    def nbonds (self):
        if hasattr (self, "bonds"):
            return len (self.bonds)
        return 0


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
