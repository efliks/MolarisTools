#-------------------------------------------------------------------------------
# . File      : MolarisInput_ToEVBTypes.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisTools.Units    import DEFAULT_EVB_LIB
from MolarisTools.Parser   import MolarisInputFile
from MolarisTools.Library  import EVBLibrary, EVBMorseAtom, EVBMorsePair


def MolarisInput_ToEVBTypes (filename, evbLibrary=DEFAULT_EVB_LIB, logging=True):
    """Write a list of EVB parameters for a given Molaris input file."""
    molarisInput = MolarisInputFile (filename, logging=logging)
    library      = EVBLibrary (evbLibrary, logging=logging)
    states       = molarisInput.types
    types        = []
    for state in states:
        for evbType in state:
            if evbType not in types:
                types.append (evbType)
    return types


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
