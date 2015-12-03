#!/usr/bin/env python

import sys, os.path
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import MolarisResidue


# . Parse the output file from Molaris and assign atomic charges and types to the EVB atoms.
# . Charges for the RS are now as in Enzymix and calculated with CHELPG.
# . Atom types are taken from Han.
nameToChargeType = {
 "PA"   :   (( 1.29 ,   "P0")  ,   ( 1.175  ,   "P+")  ,   ( 1.06 ,   "P0"))  ,#
 "O1A"  :   ((-0.85 ,   "O-")  ,   (-0.77   ,   "O-")  ,   (-0.69 ,   "O-"))  ,#
 "O2A"  :   ((-0.85 ,   "O-")  ,   (-0.77   ,   "O-")  ,   (-0.69 ,   "O-"))  ,#
 "O3A"  :   ((-0.64 ,   "Op")  ,   (-0.80   ,   "Op")  ,   (-0.96 ,   "O-"))  ,#
\
 "PB"   :   (( 1.25 ,   "P0")  ,   ( 1.19   ,   "P0")  ,   ( 1.13 ,   "P0"))  ,#
 "O1B"  :   ((-0.82 ,   "O-")  ,   (-0.89   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O2B"  :   ((-0.82 ,   "O-")  ,   (-0.89   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O3B"  :   ((-0.51 ,   "Op")  ,   (-0.505  ,   "Op")  ,   (-0.50 ,   "Op"))  ,#
\
 "PG"   :   (( 1.11 ,   "P0")  ,   ( 1.12   ,   "P0")  ,   ( 1.13 ,   "P0"))  ,#
 "O1G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O2G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O3G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
\
 "O5'"  :   ((-0.52 ,   "O0")  ,   (-0.50   ,   "O0")  ,   (-0.48 ,   "O0"))  ,#
\
 "O3'"  :   ((-0.84 ,   "O3")  ,   (-0.705  ,   "Oq")  ,   (-0.57 ,   "Oq"))  ,#
 "C3'"  :   ((-0.13 ,   "C0")  ,   ( 0.115  ,   "C0")  ,   ( 0.36 ,   "C0"))  ,#
 "H3'"  :   (( 0.14 ,   "H0")  ,   ( 0.13   ,   "H0")  ,   ( 0.12 ,   "H0"))  ,#
 "C2'"  :   ((-0.31 ,   "C0")  ,   (-0.26   ,   "C0")  ,   (-0.21 ,   "C0"))  ,#
 "H2'1" :   (( 0.07 ,   "H0")  ,   ( 0.06   ,   "H0")  ,   ( 0.05 ,   "H0"))  ,#
 "H2'2" :   (( 0.07 ,   "H0")  ,   ( 0.06   ,   "H0")  ,   ( 0.05 ,   "H0"))  ,#
}

resonanceForms    = (2, 3) # (1, 2)
verbose           = False
tab               = 2

atomsNucleotide   = ("O5'", "PA", "O1A", "O2A", "O3A", "PB", "O1B", "O2B", "O3B", "PG", "O1G", "O2G", "O3G", )
atomsPrimer       = ("O3'", "C3'", "H3'", "C2'", "H2'1", "H2'2", )
nucleotide        = MolarisResidue ("NUX" , 1 , atomsNucleotide , "determine_atoms.out")
primer            = MolarisResidue ("PRX" , 2 , atomsPrimer     , "determine_atoms.out")

nucleotide.WriteEVBAtoms (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)
primer.WriteEVBAtoms     (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)

nucleotide.WriteEVBBonds (nameToChargeType, resonanceForms, tab=tab)
primer.WriteEVBBonds     (nameToChargeType, resonanceForms, tab=tab)
