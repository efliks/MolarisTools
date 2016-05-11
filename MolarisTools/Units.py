#-------------------------------------------------------------------------------
# . File      : Units.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# http://users.mccammon.ucsd.edu/~dzhang/energy-unit-conv-table.html

HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM = 627.509469 / 0.529177       # 1185.8215096272136
HARTREE_TO_KCAL_MOL               = 627.509469
GRADIENT_TO_FORCE                 =  -1.
EV_TO_KCAL_MOL                    =  23.0609

BOHR_TO_ANGSTROM                  =   0.529177
ANGSTROM_TO_BOHR                  =   1. / BOHR_TO_ANGSTROM


atomicNumberToSymbol = {
    1   :  "H"  ,
    3   :  "Li" ,
    6   :  "C"  ,
    7   :  "N"  ,
    8   :  "O"  ,
    9   :  "F"  ,
    11  :  "NA" ,
    12  :  "MG" ,
    15  :  "P"  ,
    16  :  "S"  ,
    17  :  "CL" ,
    35  :  "BR" , }

# . Dictionary comprehension does not work in Python 2.6
symbolToAtomicNumber = {}
for number, symbol in atomicNumberToSymbol.iteritems ():
    symbolToAtomicNumber[symbol] = number
