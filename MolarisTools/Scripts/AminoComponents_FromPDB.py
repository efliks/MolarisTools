#-------------------------------------------------------------------------------
# . File      : AminoComponents_FromPDB.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import math

from MolarisTools.Units    import typicalBonds
from MolarisTools.Parser   import PDBFile
from MolarisTools.Library  import AminoComponent, AminoGroup, AminoAtom


_DEFAULT_FORCE          = 5.
_DEFAULT_TOLERANCE      = 1.6
_DEFAULT_TOLERANCE_LOW  = 0.85

# . Generally typical bond lengths seem to work better if slightly stretched
_SCALE_TOLERANCE        = 1.2


def GetBondTolerance (atoma, atomb, defaultTolerance, labels=("CL", "BR", ), logging=True):
    """Get a typical bond length between two atoms."""
    templ = atoma.label.upper ()[:2]
    if templ in labels:
        labela = templ
    else:
        labela = templ[:1]
    templ = atomb.label.upper ()[:2]
    if templ in labels:
        labelb = templ
    else:
        labelb = templ[:1]
    key = (labela, labelb)
    if typicalBonds.has_key (key):
        tolerance = typicalBonds[key]
    else:
        yek = (labelb, labela)
        if typicalBonds.has_key (yek):
            tolerance = typicalBonds[yek]
        else:
            tolerance = defaultTolerance / _SCALE_TOLERANCE
            if logging:
                print ("# . Warning: Using default tolerance for atom pair (%s, %s)" % (atoma.label, atomb.label))
    return (tolerance * _SCALE_TOLERANCE)


def BondsFromDistances (atoms, tolerance=_DEFAULT_TOLERANCE, toleranceLow=_DEFAULT_TOLERANCE_LOW, logging=True):
    """Generate a list of bonds based on distances between atoms."""
    # . Construct a distance matrix
    rows   = []
    for i, (label, serial, x, y, z) in enumerate (atoms):
        columns = []
        for j, (otherLabel, otherSerial, otherX, otherY, otherZ) in enumerate (atoms):
            if j < i:
                distance = math.sqrt ((x - otherX) ** 2 + (y - otherY) ** 2 + (z - otherZ) ** 2)
                columns.append (distance)
            elif j == i:
                columns.append (0.)
            else:
                break
        rows.append (columns)
    if logging:
        print ("# . Distance matrix")
        line = " " * 8
        for atom in atoms:
            line = "%s%6d" % (line, atom.serial)
        print line
        line = " " * 8
        for atom in atoms:
            line = "%s%6s" % (line, atom.label)
        print line
        for i, (atom, row) in enumerate (zip (atoms, rows)):
            line = "%8s" % ("%3d %4s" % (atom.serial, atom.label))
            for j, column in enumerate (row):
                mark = ""
                if i != j:
                    if column <= GetBondTolerance (atoms[i], atoms[j], tolerance, logging=False):
                        if column <= toleranceLow:
                            mark = "!"
                        else:
                            mark = "*"
                line = "%s%6s" % (line, "%s%.2f" % (mark, column))
            print line
    # . Search the distance matrix for bonds
    bonds = []
    for i, row in enumerate (rows):
        for j, column in enumerate (row):
            if i != j:
                if column <= GetBondTolerance (atoms[i], atoms[j], tolerance, logging=logging):
                    if column <= toleranceLow:
                        if logging:
                            print ("# . Warning: Atoms (%s, %s) are too close to each other" % (atoma.label, atomb.label))
                    (atoma, atomb) = (atoms[i], atoms[j])
                    pair  = ((i, atoma.serial, atoma.label), (j, atomb.serial, atomb.label))
                    bonds.append (pair)
    if logging:
        nbonds = len (bonds)
        print ("# . Found %d bonds" % nbonds)
        for i, ((indexa, seriala, labela), (indexb, serialb, labelb)) in enumerate (bonds, 1):
            print ("%3d      %4s %3d    %4s %3d" % (i, labela, seriala, labelb, serialb))
    return bonds


def AminoComponents_FromPDB (filename, tolerance=_DEFAULT_TOLERANCE, toleranceLow=_DEFAULT_TOLERANCE_LOW, logging=True, verbose=True):
    """Automatically generate components based on a PDB file."""
    pdb        = PDBFile (filename)
    components = []
    for residue in pdb.residues:
        if logging:
            print ("*** Building component: %s %s %d ***" % (residue.chain, residue.label, residue.serial))
        # . Generate unique atom labels
        uniqueLabels = []
        i = 1
        for atom in residue.atoms:
            label = atom.label
            while label in uniqueLabels:
                label = "%s%d" % (label, i)
                i += 1
            uniqueLabels.append (label)
        # . Generate atoms
        aminoAtoms = []
        for atom, uniqueLabel in zip (residue.atoms, uniqueLabels):
            atype = "%1s0" % atom.label[0]
            aminoAtom = AminoAtom (
                atomLabel  = uniqueLabel ,
                atomType   = atype       ,
                atomCharge = 0.          ,
                )
            aminoAtoms.append (aminoAtom)
        # . Try to read bonds from the PDB file
        bonds = residue.GetBonds ()
        if not bonds:
            # . Generate bonds
            bonds = BondsFromDistances (residue.atoms, logging=logging, tolerance=tolerance, toleranceLow=toleranceLow)
        aminoBonds  = []
        for (i, seriala, labela), (j, serialb, labelb) in bonds:
            pair    = (uniqueLabels[i], uniqueLabels[j])
            aminoBonds.append (pair)
        # . Generate groups
        natoms     = len (residue.atoms)
        aminoGroup = AminoGroup (
            natoms       =  natoms           ,
            centralAtom  =  uniqueLabels[0]  ,
            labels       =  uniqueLabels     ,
            radius       =  3.               ,
            symbol       =  "A"              ,
            )
        aminoGroups = [aminoGroup, ]
        # . Finally, generate a component
        component = AminoComponent (
            serial  =   residue.serial  ,
            name    =   residue.label   ,
            atoms   =   aminoAtoms      ,
            bonds   =   aminoBonds      ,
            groups  =   aminoGroups     ,
            connect =   ("", "")        ,
            logging =   logging         ,
            verbose =   verbose         ,
            title   =   "Automatically generated from %s %s %d" % (residue.chain, residue.label, residue.serial) ,
            )
        components.append (component)
    return components


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
