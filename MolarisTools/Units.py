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

atomicNumberToSymbol = {1 : "H", 6 : "C", 7 : "N", 8 : "O", 9 : "F", 17 : "CL", 15 : "P", 12 : "MG", 35 : "BR"}

# . Dictionary comprehension does not work for Python 2.6
symbolToAtomicNumber = {}
for number, symbol in atomicNumberToSymbol.iteritems ():
    symbolToAtomicNumber[symbol] = number
