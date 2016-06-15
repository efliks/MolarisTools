#-------------------------------------------------------------------------------
# . File      : Scripts.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisOutputFile      import MolarisOutputFile
from AminoLibrary           import AminoLibrary, AminoComponent, AminoGroup, AminoAtom
from PDBFile                import PDBFile

import os, math


_DEFAULT_LIBRARY_FILE   = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "amino98_custom_small.lib")
_DEFAULT_FORCE          = 5.

_DEFAULT_TOLERANCE      = 1.6
_DEFAULT_TOLERANCE_LOW  = 0.85


def GenerateEVBList (fileLibrary=_DEFAULT_LIBRARY_FILE, fileMolarisOutput="determine_atoms.out", selectGroups={}, ntab=2, exceptions=("MG", "CL", )):
    """Generate a list of EVB atoms and bonds based on a Molaris output file."""
    library = AminoLibrary (fileLibrary, logging=False)
    
    mof        = MolarisOutputFile (fileMolarisOutput)
    tabs       = (" " * 4) * ntab
    natoms     = 0
    nbonds     = 0
    evbSerials = []
    for residue in mof.residues:
        libResidue    = library[residue.label]
        print ("%s# . %s%d" % (tabs, residue.label, residue.serial))
    
        # . Print EVB atoms
        #          evb_atm         6   -0.3000    C0   -0.3000    C0      #    -0.3000  CT    C5' 
        #          evb_atm         5   -0.4500    O0   -0.4500    O0      #    -0.4500  O4    O5' 
        #   (...)
        for atom in residue.atoms:
            if atom.atype in exceptions:
                evbType = atom.atype
            else:
                evbType = "%1s0" % atom.atype[0]
            groupLabel = "?"
            for group in libResidue.groups:
                if atom.label in group.labels:
                    groupLabel = group.symbol
                    break
            # . If selectGroups is empty, select all atoms from all residues.
            # . Otherwise, select only atoms belonging to specific groups in residues.
            #
            # . selectGroups has the following format:
            #       "residueLabel" : ("groupLabel1", "groupLabel2", )
            includeAtom = True
            if selectGroups != {}:
                includeAtom = False
                if selectGroups.has_key (residue.label):
                    groupLabels = selectGroups[residue.label]
                    if groupLabel in groupLabels:
                        includeAtom = True
            if includeAtom:
                print ("%sevb_atm    %2d    %6.3f    %2s        %6.3f    %2s    #  %6.3f  %4s  %4s    %1s" % (tabs, atom.serial, atom.charge, evbType, atom.charge, evbType, atom.charge, atom.atype, atom.label, groupLabel))
                evbSerials.append (atom.serial)
                natoms += 1
    
        # . Print EVB bonds
        #          evb_bnd   0         5     6   # O5'     C5'   
        #          evb_bnd   0         7     6   # H5'1    C5'   
        #   (...)
        pairs  = []
        labels = []
        for atom in residue.atoms:
            # . Check if the atom is an EVB atom
            if atom.serial in evbSerials:
                for (serial, label) in atom.bonds:
                    # . Connect only to another EVB atom
                    if serial in evbSerials:
                        pair    = (atom.serial, serial     )
                        pairRev = (serial     , atom.serial)
                        if not ((pair in pairs) or (pairRev in pairs)):
                            pairs.append (pair)
                            pairLabels = (atom.label, label)
                            labels.append (pairLabels)
        for (seriala, serialb), (labela, labelb) in zip (pairs, labels):
            print ("%sevb_bnd   0   %3d   %3d    #  %3s  %3s" % (tabs, seriala, serialb, labela, labelb))
            nbonds += 1
    # . Print summary
    print ("%s# EVB atoms: %d, EVB bonds: %d" % (tabs, natoms, nbonds))

    # . Print constrained atoms
    # constraint_post     1     5.0   5.0   5.0     5.138     5.342    14.879    add_to_qm_force    0   # PG
    #   (...)
    tabs = (" " * 4) * (ntab + 1)
    for residue in mof.residues:
        print ("%s# . %s%d" % (tabs, residue.label, residue.serial))
        for atom in residue.atoms:
            if atom.serial in evbSerials:
                print ("%s# constraint_post  %4d    %4.1f  %4.1f  %4.1f  %8.3f  %8.3f  %8.3f    add_to_qm_force    0   # %s" % (tabs, atom.serial, _DEFAULT_FORCE, _DEFAULT_FORCE, _DEFAULT_FORCE, atom.x, atom.y, atom.z, atom.label))


def DetermineBAT (fileLibrary=_DEFAULT_LIBRARY_FILE, fileMolarisOutput="determine_atoms.out", residueLabels=(), ):
    """Determine bonds, angles and torsions to analyze their statistical distributions."""

    # . Load the log file
    mof     = MolarisOutputFile (fileMolarisOutput)
    # . Load the library
    library = AminoLibrary (fileLibrary, logging=False)

    if mof.nresidues > 0:
        for residue in mof.residues:
            include = True
            if len (residueLabels) > 0:
                if residue.label not in residueLabels:
                    include = False
            if include:
                component = library[residue.label]
                component.GenerateAngles (quiet=True)
                component.GenerateTorsions (quiet=True)

                # . Write bond serials
                for (bonda, bondb) in component.bonds:
                    serials = []
                    for bond in (bonda, bondb):
                        for atom in residue.atoms:
                            if atom.label == bond: break
                        serials.append (atom.serial)
                    print ("%4d    %4d    # %4s    %4s" % (serials[0], serials[1], bonda, bondb))

                # . Write angle serials
                for (anglea, angleb, anglec) in component.angles:
                    serials = []
                    for angle in (anglea, angleb, anglec):
                        for atom in residue.atoms:
                            if atom.label == angle: break
                        serials.append (atom.serial)
                    print ("%4d    %4d    %4d    # %4s    %4s    %4s" % (serials[0], serials[1], serials[2], anglea, angleb, anglec))

                # . Write torsion serials
                for (torsiona, torsionb, torsionc, torsiond) in component.torsions:
                    serials = []
                    for torsion in (torsiona, torsionb, torsionc, torsiond):
                        for atom in residue.atoms:
                            if atom.label == torsion: break
                        serials.append (atom.serial)
                    print ("%4d    %4d    %4d    %4d    # %4s    %4s    %4s    %4s" % (serials[0], serials[1], serials[2], serials[3], torsiona, torsionb, torsionc, torsiond))


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
        for atom, row in zip (atoms, rows):
            line = "%8s" % ("%3d %4s" % (atom.serial, atom.label))
            for column in row:
                line = "%s%6.2f" % (line, column)
            print line
    # . Search the distance matrix for bonds
    # . TODO
    # . Have a table of common bond lengths and replace tolerance by a percentage of
    # . how much a bond can differ from its tabled value.
    bonds = []
    for i, row in enumerate (rows):
        for j, column in enumerate (row):
            if i != j:
                if column <= tolerance:
                    atoma = atoms[i]
                    atomb = atoms[j]
                    if column <= toleranceLow:
                        if logging:
                            print ("# . Warning: Atoms %s%d and %s%d are too close to each other." % (atoma.label, atoma.serial, atomb.label, atomb.serial))
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
        # . Generate bonds
        bonds       = BondsFromDistances (residue.atoms, logging=logging)
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
if __name__ == "__main__": pass
