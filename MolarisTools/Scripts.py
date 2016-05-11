#-------------------------------------------------------------------------------
# . File      : Scripts.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisOutputFile      import MolarisOutputFile
from AminoLibrary           import AminoLibrary

import os


_DEFAULT_LIBRARY_FILE = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "amino98_custom_small.lib")

def GenerateEVBList (fileLibrary=_DEFAULT_LIBRARY_FILE, fileMolarisOutput="determine_atoms.out", ntab=2):
    """Generate a list of EVB atoms and bonds based on a Molaris output file."""
    library = AminoLibrary (fileLibrary, logging=False)
    
    mof     = MolarisOutputFile (fileMolarisOutput)
    tabs    = (" " * 4) * ntab
    natoms  = 0
    nbonds  = 0
    for residue in mof.residues:
        libResidue = library[residue.label]
        print ("%s# . %s%d" % (tabs, residue.label, residue.serial))
    
        # . Print EVB atoms
        #          evb_atm         6   -0.3000    C0   -0.3000    C0      #    -0.3000  CT    C5' 
        #          evb_atm         5   -0.4500    O0   -0.4500    O0      #    -0.4500  O4    O5' 
        #   (...)
        for atom in residue.atoms:
            evbType = "%1s0" % atom.atype[0]
            groupLabel = "?"
            for group in libResidue.groups:
                if atom.label in group.labels:
                    groupLabel = group.symbol
                    break
            print ("%sevb_atm    %2d    %6.3f    %2s        %6.3f    %2s    #  %6.3f  %2s  %2s    %2s" % (tabs, atom.serial, atom.charge, evbType, atom.charge, evbType, atom.charge, atom.atype, atom.label, groupLabel))
            natoms += 1
    
        # . Print EVB bonds
        #          evb_bnd   0         5     6   # O5'     C5'   
        #          evb_bnd   0         7     6   # H5'1    C5'   
        #   (...)
        pairs  = []
        labels = []
        for atom in residue.atoms:
            for (serial, label) in atom.bonds:
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


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
