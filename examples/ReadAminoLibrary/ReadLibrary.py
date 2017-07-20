#-------------------------------------------------------------------------------
# . File      : ReadAmino.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisTools.Library  import AminoLibrary


library = AminoLibrary (logging=True, verbose=True, unique=True, filename="amino98_tiny.lib")

component = library["XPR"]
component.Write (filename="component.lib", terminate=True)
