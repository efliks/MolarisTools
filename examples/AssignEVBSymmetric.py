#!/usr/bin/env python

import sys, os.path
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import MolarisResidue


# . Parse the output file from Molaris and assign atomic charges and types to the EVB atoms.
# . Charges for the RS are now as in Enzymix and calculated with CHELPG.
# . Atom types are taken from Han.
# . TS1, INT and TS2 share the same atom types.
#                  RS                   TS1                 INT                 TS2                 PS
nameToChargeType = {
# NUS phosphate
"PA"    :  (( 1.24 , "P0" ) ,    ( 1.30   , "P+" ) ,    ( 1.30  , "P+" ) ,    ( 1.30   , "P+" ) ,    ( 1.24 , "P0" )) ,
"O1A"   :  ((-0.71 , "O-" ) ,    (-0.74   , "O-" ) ,    (-0.74  , "O-" ) ,    (-0.74   , "O-" ) ,    (-0.71 , "O-" )) ,
"O2A"   :  ((-0.71 , "O-" ) ,    (-0.74   , "O-" ) ,    (-0.74  , "O-" ) ,    (-0.74   , "O-" ) ,    (-0.71 , "O-" )) ,
# NUS side-group
"O5'"   :  ((-0.41 , "Oq" ) ,    (-0.41   , "Oq" ) ,    (-0.41  , "Oq" ) ,    (-0.41   , "Oq" ) ,    (-0.41 , "Oq" )) ,
"C5'"   :  ((-0.24 , "C0" ) ,    (-0.24   , "C0" ) ,    (-0.24  , "C0" ) ,    (-0.24   , "C0" ) ,    (-0.24 , "C0" )) ,
"H5'1"  :  (( 0.08 , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08  , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08 , "H0" )) ,
"H5'2"  :  (( 0.08 , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08  , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08 , "H0" )) ,
"H5'3"  :  (( 0.08 , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08  , "H0" ) ,    ( 0.08   , "H0" ) ,    ( 0.08 , "H0" )) ,
# NUS leaving group
"O3A"   :  ((-0.41 , "Oq" ) ,    (-0.5325 , "Oq" ) ,    (-0.655 , "Oq" ) ,    (-0.7775 , "Oq" ) ,    (-0.90 , "O3" )) ,
"C6'"   :  ((-0.24 , "C0" ) ,    (-0.2425 , "C0" ) ,    (-0.245 , "C0" ) ,    (-0.2475 , "C0" ) ,    (-0.25 , "C0" )) ,
"H6'1"  :  (( 0.08 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.05 , "H0" )) ,
"H6'2"  :  (( 0.08 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.05 , "H0" )) ,
"H6'3"  :  (( 0.08 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.05 , "H0" )) ,
# PRS attacking group
"O3'"   :  ((-0.90 , "O3" ) ,    (-0.7775 , "Oq" ) ,    (-0.655 , "Oq" ) ,    (-0.5325 , "Oq" ) ,    (-0.41 , "Oq" )) ,
"C3'"   :  ((-0.25 , "C0" ) ,    (-0.2475 , "C0" ) ,    (-0.245 , "C0" ) ,    (-0.2425 , "C0" ) ,    (-0.24 , "C0" )) ,
"H3'1"  :  (( 0.05 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.08 , "H0" )) ,
"H3'2"  :  (( 0.05 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.08 , "H0" )) ,
"H3'3"  :  (( 0.05 , "H0" ) ,    ( 0.0575 , "H0" ) ,    ( 0.065 , "H0" ) ,    ( 0.0725 , "H0" ) ,    ( 0.08 , "H0" )) ,
}
#                   H6'2
#                    |
#              H6'1-C6'-H6'3
#      H5'1         /
#        |         /
#  H5'2-C5'      O3A
#       / \       |  O2A(-)
#    H5'3  \      | /
#           O5'---PA
#                 | \\
#                 |  O1A
#             (-)O3'
#                /
#               /
#         H3'3-C3'-H3'1
#              |
#             H3'2

stateToNumber = { "RS"  : 1 ,
                  "TS1" : 2 ,
                  "INT" : 3 ,
                  "TS2" : 4 ,
                  "PS"  : 5 , }
resonanceForms    = (stateToNumber["RS"], stateToNumber["TS1"])
verbose           = False
tab               = 2

atomsNucleotide   = ("PA", "O1A", "O2A", "O3A", "O5'", "C5'", \
    "H5'1", "H5'2", "H5'3", "C6'", "H6'1", "H6'2", "H6'3")
atomsPrimer       = ("O3'", "C3'" , "H3'1", "H3'2", "H3'3", )
nucleotide        = MolarisResidue ("NUS" , 1 , atomsNucleotide , "./determine_atoms.out")
primer            = MolarisResidue ("PRS" , 2 , atomsPrimer     , "./determine_atoms.out")

nucleotide.WriteEVBAtoms (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)
primer.WriteEVBAtoms     (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)

nucleotide.WriteEVBBonds (nameToChargeType, resonanceForms, tab=tab)
primer.WriteEVBBonds     (nameToChargeType, resonanceForms, tab=tab)
