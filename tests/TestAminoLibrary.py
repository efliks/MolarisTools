#-------------------------------------------------------------------------------
# . File      : TestAminoLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
# . Test 01   : Operations on amino-library
#-------------------------------------------------------------------------------
import unittest, sys, os

from MolarisTools.Library  import AminoLibrary


class TestAminoLibrary (unittest.TestCase):
    def test_Read (self):
        library  = AminoLibrary (filename=os.path.join ("..", "data", "amino98_custom_small.lib"), logging=True, verbose=False)
        self.assertEqual (len (library), 153)


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"):
    unittest.main ()
