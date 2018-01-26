#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : CallTeraChem.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import sys, os

pathMolarisTools = os.path.join (os.environ["HOME"], "devel", "MolarisTools")
sys.path.append (pathMolarisTools)

from MolarisTools.QMMM  import QMCallerTeraChem


callerTeraChem = QMCallerTeraChem (
        qmmm            =  False           ,
        archive         =  False           ,
        restart         =  True            ,
        charge          =   -1             ,
        multiplicity    =    1             ,
        method          =  "B3LYP/6-31G*"  ,
        pathTeraChem    = os.path.join (os.environ["HOME"], "TeraChem", "bin", "terachem") ,
            )
callerTeraChem.Run ()

# chargeScheme = "MerzKollman"