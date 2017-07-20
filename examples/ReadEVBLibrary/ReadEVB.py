#-------------------------------------------------------------------------------
# . File      : ReadEVB.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisTools.Library  import EVBLibrary


lib = EVBLibrary (filename="evb_poll_clean.lib")

# . Reduce the library to selected atom types.
types = ("L-", "L0", 
         "B-", "B0", 
         "C0", "C+", "H0", )
lib.PurgeTypes (types)

lib.WriteLibrary (digits=2)
