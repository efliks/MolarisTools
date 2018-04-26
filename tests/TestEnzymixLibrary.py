#-------------------------------------------------------------------------------
# . File      : TestEnzymixLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
# . Test 03   : Operations on ENZYMIX library
#-------------------------------------------------------------------------------
import unittest, sys, os

from MolarisTools.Library  import ParametersLibrary


class TestEnzymixLibrary (unittest.TestCase):
    def test_Read (self):
        library  = ParametersLibrary (filename=os.path.join ("..", "data", "parm.lib"), logging=True)
        self.assertEqual (len (library), 815)


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"):
    unittest.main ()
