#-------------------------------------------------------------------------------
# . File      : Units.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os

DEFAULT_AMINO_LIB  = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "amino98_custom_small.lib")
DEFAULT_PARM_LIB   = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "parm.lib")
DEFAULT_EVB_LIB    = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "evb_poll_clean.lib")


# . http://users.mccammon.ucsd.edu/~dzhang/energy-unit-conv-table.html
BOHR_TO_ANGSTROM                  =   0.529177
ANGSTROM_TO_BOHR                  =   1. / BOHR_TO_ANGSTROM
HARTREE_TO_KCAL_MOL               = 627.509469
HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM = HARTREE_TO_KCAL_MOL / BOHR_TO_ANGSTROM       # 1185.8215096272136
EV_TO_KCAL_MOL                    =  23.0609
GRADIENT_TO_FORCE                 =  -1.


# . Typical bond lengths
# . http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
typicalBonds = {
    ("H"   ,  "H" )    :   0.74  ,
    ("H"   ,  "F" )    :   0.92  ,
    ("H"   ,  "O" )    :   0.96  ,
    ("H"   ,  "N" )    :   1.01  ,
    ("H"   ,  "C" )    :   1.09  ,
    ("H"   ,  "CL")    :   1.27  ,
    ("C"   ,  "F" )    :   1.35  ,
    ("H"   ,  "BR")    :   1.41  ,
    ("F"   ,  "F" )    :   1.42  ,
    ("C"   ,  "O" )    :   1.43  ,
    ("N"   ,  "N" )    :   1.45  ,
    ("C"   ,  "N" )    :   1.47  ,
    ("O"   ,  "O" )    :   1.48  ,
    ("C"   ,  "C" )    :   1.54  ,
    ("P"   ,  "F" )    :   1.54  ,
    ("H"   ,  "I" )    :   1.61  ,
    ("P"   ,  "O" )    :   1.63  ,
    ("C"   ,  "CL")    :   1.77  ,
    ("C"   ,  "S" )    :   1.82  ,
    ("P"   ,  "C" )    :   1.84  ,
    ("P"   ,  "S" )    :   1.86  ,
    ("C"   ,  "BR")    :   1.94  ,
    ("CL"  ,  "CL")    :   1.99  ,
    ("P"   ,  "CL")    :   2.03  ,
    ("C"   ,  "I" )    :   2.14  ,
    ("P"   ,  "P" )    :   2.21  ,
    ("BR"  ,  "BR")    :   2.28  ,
    ("I"   ,  "I" )    :   2.67  ,
}


atomicNumberToSymbol = {
   1   :  "H" ,     #  Hydrogen
   2   :  "HE" ,    #  Helium
   3   :  "LI" ,    #  Lithium
   4   :  "BE" ,    #  Beryllium
   5   :  "B" ,     #  Boron
   6   :  "C" ,     #  Carbon
   7   :  "N" ,     #  Nitrogen
   8   :  "O" ,     #  Oxygen
   9   :  "F" ,     #  Fluorine
  10   :  "NE" ,    #  Neon
  11   :  "NA" ,    #  Sodium
  12   :  "MG" ,    #  Magnesium
  13   :  "AL" ,    #  Aluminum
  14   :  "SI" ,    #  Silicon
  15   :  "P" ,     #  Phosphorus
  16   :  "S" ,     #  Sulfur
  17   :  "CL" ,    #  Chlorine
  18   :  "AR" ,    #  Argon
  19   :  "K" ,     #  Potassium
  20   :  "CA" ,    #  Calcium
  21   :  "SC" ,    #  Scandium
  22   :  "TI" ,    #  Titanium
  23   :  "V" ,     #  Vanadium
  24   :  "CR" ,    #  Chromium
  25   :  "MN" ,    #  Manganese
  26   :  "FE" ,    #  Iron
  27   :  "CO" ,    #  Cobalt
  28   :  "NI" ,    #  Nickel
  29   :  "CU" ,    #  Copper
  30   :  "ZN" ,    #  Zinc
  31   :  "GA" ,    #  Gallium
  32   :  "GE" ,    #  Germanium
  33   :  "AS" ,    #  Arsenic
  34   :  "SE" ,    #  Selenium
  35   :  "BR" ,    #  Bromine
  36   :  "KR" ,    #  Krypton
  37   :  "RB" ,    #  Rubidium
  38   :  "SR" ,    #  Strontium
  39   :  "Y" ,     #  Yttrium
  40   :  "ZR" ,    #  Zirconium
  41   :  "NB" ,    #  Niobium
  42   :  "MO" ,    #  Molybdenum
  43   :  "TC" ,    #  Technetium
  44   :  "RU" ,    #  Ruthenium
  45   :  "RH" ,    #  Rhodium
  46   :  "PD" ,    #  Palladium
  47   :  "AG" ,    #  Silver
  48   :  "CD" ,    #  Cadmium
  49   :  "IN" ,    #  Indium
  50   :  "SN" ,    #  Tin
  51   :  "SB" ,    #  Antimony
  52   :  "TE" ,    #  Tellurium
  53   :  "I" ,     #  Iodine
  54   :  "XE" ,    #  Xenon
  55   :  "CS" ,    #  Cesium
  56   :  "BA" ,    #  Barium
  57   :  "LA" ,    #  Lanthanum
  58   :  "CE" ,    #  Cerium
  59   :  "PR" ,    #  Praseodymium
  60   :  "ND" ,    #  Neodymium
  61   :  "PM" ,    #  Promethium
  62   :  "SM" ,    #  Samarium
  63   :  "EU" ,    #  Europium
  64   :  "GD" ,    #  Gadolinium
  65   :  "TB" ,    #  Terbium
  66   :  "DY" ,    #  Dysprosium
  67   :  "HO" ,    #  Holmium
  68   :  "ER" ,    #  Erbium
  69   :  "TM" ,    #  Thulium
  70   :  "YB" ,    #  Ytterbium
  71   :  "LU" ,    #  Lutetium
  72   :  "HF" ,    #  Hafnium
  73   :  "TA" ,    #  Tantalum
  74   :  "W" ,     #  Tungsten
  75   :  "RE" ,    #  Rhenium
  76   :  "OS" ,    #  Osmium
  77   :  "IR" ,    #  Iridium
  78   :  "PT" ,    #  Platinum
  79   :  "AU" ,    #  Gold
  80   :  "HG" ,    #  Mercury
  81   :  "TL" ,    #  Thallium
  82   :  "PB" ,    #  Lead
  83   :  "BI" ,    #  Bismuth
  84   :  "PO" ,    #  Polonium
  85   :  "AT" ,    #  Astatine
  86   :  "RN" ,    #  Radon
  87   :  "FR" ,    #  Francium
  88   :  "RA" ,    #  Radium
  89   :  "AC" ,    #  Actinium
  90   :  "TH" ,    #  Thorium
  91   :  "PA" ,    #  Protactinium
  92   :  "U" ,     #  Uranium
  93   :  "NP" ,    #  Neptunium
  94   :  "PU" ,    #  Plutonium
  95   :  "AM" ,    #  Americium
  96   :  "CM" ,    #  Curium
  97   :  "BK" ,    #  Berkelium
  98   :  "CF" ,    #  Californium
  99   :  "ES" ,    #  Einsteinium
 100   :  "FM" ,    #  Fermium
 101   :  "MD" ,    #  Mendelevium
 102   :  "NO" ,    #  Nobelium
 103   :  "LR" ,    #  Lawrencium
 104   :  "RF" ,    #  Rutherfordium
 105   :  "DB" ,    #  Dubnium
 106   :  "SG" ,    #  Seaborgium
 107   :  "BH" ,    #  Bohrium
 108   :  "HS" ,    #  Hassium
 109   :  "MT" ,    #  Meitnerium
 110   :  "DS" ,    #  Darmstadtium
 111   :  "RG" ,    #  Roentgenium
 112   :  "UUB" ,   #  Ununbium
 113   :  "UUT" ,   #  Ununtrium
 114   :  "UUQ" ,   #  Ununquadium
 115   :  "UUP" ,   #  Ununpentium
 116   :  "UUH" ,   #  Ununhexium
 117   :  "UUS" ,   #  Ununseptium
 118   :  "UUO" ,   #  Ununoctium
    }

symbolToAtomicNumber = {}
for (number, symbol) in atomicNumberToSymbol.iteritems ():
    symbolToAtomicNumber[symbol] = number
