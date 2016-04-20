#-------------------------------------------------------------------------------
# . File      : Units.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# http://users.mccammon.ucsd.edu/~dzhang/energy-unit-conv-table.html

HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM = 627.509469 / 0.529177
HARTREE_TO_KCAL_MOL               = 627.509469
GRADIENT_TO_FORCE                 =  -1.
EV_TO_KCAL_MOL                    =  23.0609

# . Taken from pDynamo
UNITS_LENGTH_ANGSTROMS_TO_BOHRS = 1.0e-10 / 5.291772083e-11
UNITS_LENGTH_BOHRS_TO_ANGSTROMS = 5.291772083e-11 / 1.0e-10
UNITS_ENERGY_HARTREES_TO_ELECTRON_VOLTS                   = 27.2113845
UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE              = 2625.5
UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE = 4.184


atomicNumberToSymbol = {
    1   :  "H"  ,
    6   :  "C"  ,
    7   :  "N"  ,
    8   :  "O"  ,
    9   :  "F"  ,
    12  :  "MG" ,
    15  :  "P"  ,
    16  :  "S"  ,
    17  :  "CL" ,
    35  :  "BR" , }

# . Dictionary comprehension does not work in Python 2.6
symbolToAtomicNumber = {}
for number, symbol in atomicNumberToSymbol.iteritems ():
    symbolToAtomicNumber[symbol] = number
