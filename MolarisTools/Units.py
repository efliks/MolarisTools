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
   1   :  "H" ,     #      Hydrogen
   2   :  "HE" ,    #        Helium
   3   :  "LI" ,    #       Lithium
   4   :  "BE" ,    #     Beryllium
   5   :  "B" ,     #         Boron
   6   :  "C" ,     #        Carbon
   7   :  "N" ,     #      Nitrogen
   8   :  "O" ,     #        Oxygen
   9   :  "F" ,     #      Fluorine
  10   :  "NE" ,    #          Neon
  11   :  "NA" ,    #        Sodium
  12   :  "MG" ,    #     Magnesium
  13   :  "AL" ,    #      Aluminum
  14   :  "SI" ,    #       Silicon
  15   :  "P" ,     #    Phosphorus
  16   :  "S" ,     #        Sulfur
  17   :  "CL" ,    #      Chlorine
  18   :  "AR" ,    #         Argon
  19   :  "K" ,     #     Potassium
  20   :  "CA" ,    #       Calcium
  21   :  "SC" ,    #      Scandium
  22   :  "TI" ,    #      Titanium
  23   :  "V" ,     #      Vanadium
  24   :  "CR" ,    #      Chromium
  25   :  "MN" ,    #     Manganese
  26   :  "FE" ,    #          Iron
  27   :  "CO" ,    #        Cobalt
  28   :  "NI" ,    #        Nickel
  29   :  "CU" ,    #        Copper
  30   :  "ZN" ,    #          Zinc
  31   :  "GA" ,    #       Gallium
  32   :  "GE" ,    #     Germanium
  33   :  "AS" ,    #       Arsenic
  34   :  "SE" ,    #      Selenium
  35   :  "BR" ,    #       Bromine
  36   :  "KR" ,    #       Krypton
  37   :  "RB" ,    #      Rubidium
  38   :  "SR" ,    #     Strontium
  39   :  "Y" ,     #       Yttrium
  40   :  "ZR" ,    #     Zirconium
  41   :  "NB" ,    #       Niobium
  42   :  "MO" ,    #    Molybdenum
  43   :  "TC" ,    #    Technetium
  44   :  "RU" ,    #     Ruthenium
  45   :  "RH" ,    #       Rhodium
  46   :  "PD" ,    #     Palladium
  47   :  "AG" ,    #        Silver
  48   :  "CD" ,    #       Cadmium
  49   :  "IN" ,    #        Indium
  50   :  "SN" ,    #           Tin
  51   :  "SB" ,    #      Antimony
  52   :  "TE" ,    #     Tellurium
  53   :  "I" ,     #        Iodine
  54   :  "XE" ,    #         Xenon
  55   :  "CS" ,    #        Cesium
  56   :  "BA" ,    #        Barium
  57   :  "LA" ,    #     Lanthanum
  58   :  "CE" ,    #        Cerium
  59   :  "PR" ,    #  Praseodymium
  60   :  "ND" ,    #     Neodymium
  61   :  "PM" ,    #    Promethium
  62   :  "SM" ,    #      Samarium
  63   :  "EU" ,    #      Europium
  64   :  "GD" ,    #    Gadolinium
  65   :  "TB" ,    #       Terbium
  66   :  "DY" ,    #    Dysprosium
  67   :  "HO" ,    #       Holmium
  68   :  "ER" ,    #        Erbium
  69   :  "TM" ,    #       Thulium
  70   :  "YB" ,    #     Ytterbium
  71   :  "LU" ,    #      Lutetium
  72   :  "HF" ,    #       Hafnium
  73   :  "TA" ,    #      Tantalum
  74   :  "W" ,     #      Tungsten
  75   :  "RE" ,    #       Rhenium
  76   :  "OS" ,    #        Osmium
  77   :  "IR" ,    #       Iridium
  78   :  "PT" ,    #      Platinum
  79   :  "AU" ,    #          Gold
  80   :  "HG" ,    #       Mercury
  81   :  "TL" ,    #      Thallium
  82   :  "PB" ,    #          Lead
  83   :  "BI" ,    #       Bismuth
  84   :  "PO" ,    #      Polonium
  85   :  "AT" ,    #      Astatine
  86   :  "RN" ,    #         Radon
  87   :  "FR" ,    #      Francium
  88   :  "RA" ,    #        Radium
  89   :  "AC" ,    #      Actinium
  90   :  "TH" ,    #       Thorium
  91   :  "PA" ,    #  Protactinium
  92   :  "U" ,     #       Uranium
  93   :  "NP" ,    #     Neptunium
  94   :  "PU" ,    #     Plutonium
  95   :  "AM" ,    #     Americium
  96   :  "CM" ,    #        Curium
  97   :  "BK" ,    #     Berkelium
  98   :  "CF" ,    #   Californium
  99   :  "ES" ,    #   Einsteinium
 100   :  "FM" ,    #       Fermium
 101   :  "MD" ,    #   Mendelevium
 102   :  "NO" ,    #      Nobelium
 103   :  "LR" ,    #    Lawrencium
 104   :  "RF" ,    # Rutherfordium
 105   :  "DB" ,    #       Dubnium
 106   :  "SG" ,    #    Seaborgium
 107   :  "BH" ,    #       Bohrium
 108   :  "HS" ,    #       Hassium
 109   :  "MT" ,    #    Meitnerium
 110   :  "DS" ,    #  Darmstadtium
 111   :  "RG" ,    #   Roentgenium
 112   :  "UUB" ,   #      Ununbium
 113   :  "UUT" ,   #     Ununtrium
 114   :  "UUQ" ,   #   Ununquadium
 115   :  "UUP" ,   #   Ununpentium
 116   :  "UUH" ,   #    Ununhexium
 117   :  "UUS" ,   #   Ununseptium
 118   :  "UUO" ,   #    Ununoctium
    }


# . Dictionary comprehension does not work in Python 2.6
symbolToAtomicNumber = {}
for number, symbol in atomicNumberToSymbol.iteritems ():
    symbolToAtomicNumber[symbol] = number
