#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : CallMopac.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os

from MolarisTools.QMMM  import QMCallerMopac


# . Modify the electronic state and method.
callerMopac = QMCallerMopac (
        archive         =  True     ,
        charge          =    0      ,
        multiplicity    =    1      ,
        method          =  "MNDOD"  ,
        qmmm            =  True     ,
        pathMopac       = os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe") ,
            )
callerMopac.Run ()
