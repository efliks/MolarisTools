#!/usr/bin/env python
#-------------------------------------------------------------------------------
# . File      : caller_mopac.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import sys, os

# . Define a path to the MolarisTools library or comment this out and 
# . place MolarisTools in the current directory
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))


from MolarisTools import QMCallerMopac

# . Modify the electronic state and method
callerMopac = QMCallerMopac (
        archive         =  True     ,
        charge          =    0      ,
        multiplicity    =    1      ,
        method          =  "MNDOD"  ,
        qmmm            =  True     ,
        pathMopac       = os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
            )
callerMopac.Run ()
