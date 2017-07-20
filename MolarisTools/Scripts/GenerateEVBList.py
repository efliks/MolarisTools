#-------------------------------------------------------------------------------
# . File      : GenerateEVBList.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisTools.Units    import DEFAULT_AMINO_LIB, DEFAULT_PARM_LIB, DEFAULT_EVB_LIB, typicalBonds
from MolarisTools.Parser   import MolarisOutputFile
from MolarisTools.Library  import AminoLibrary

_DEFAULT_FORCE  = 5.


def GenerateEVBList (fileLibrary=DEFAULT_AMINO_LIB, fileMolarisOutput="determine_atoms.out", selectGroups={}, ntab=2, exceptions=("MG", "CL", "BR", "DE", ), overwriteCharges=[], constrainForce=_DEFAULT_FORCE, constrainAll=False, overwriteConstraints=[]):
    """Generate a list of EVB atoms and bonds based on a Molaris output file."""
    library    = AminoLibrary (fileLibrary, logging=False)
    
    mof        = MolarisOutputFile (fileMolarisOutput)
    tabs       = (" " * 4) * ntab
    natoms     = 0
    nbonds     = 0
    evbSerials = []
    if overwriteCharges:
        charges     = iter (overwriteCharges)
    if overwriteConstraints:
        # . Do not use atomic coordinates from the PDB file, use a predefined list of positional constraints
        constraints = iter (overwriteConstraints)

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
            charge = atom.charge
            if includeAtom:
                if overwriteCharges != []:
                    charge = next (charges)
                print ("%sevb_atm    %2d    %6.2f    %2s        %6.2f    %2s    #  %6.2f  %4s  %4s    %1s" % (tabs, atom.serial, charge, evbType, charge, evbType, atom.charge, atom.atype, atom.label, groupLabel))
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
                if not constrainAll:
                    if atom.label[0] == "H":
                        continue
                (cx, cy, cz) = (atom.x, atom.y, atom.z)
                if overwriteConstraints != []:
                    (cx, cy, cz) = next (constraints)
                print ("%s# constraint_post  %4d    %4.1f  %4.1f  %4.1f  %8.3f  %8.3f  %8.3f   # %s" % (tabs, atom.serial, constrainForce, constrainForce, constrainForce, cx, cy, cz, atom.label))


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
