#-------------------------------------------------------------------------------
# . File      : TestEVBLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
# . Test 02   : Operations on EVB library
#-------------------------------------------------------------------------------
import unittest, sys, os

from MolarisTools.Library  import EVBLibrary


class TestEVBLibrary (unittest.TestCase):
    def test_Read (self):
        library  = EVBLibrary (filename=os.path.join ("..", "data", "evb_poll_clean.lib"), logging=True)
        self.assertEqual (len (library), 269)


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"):
    unittest.main ()
