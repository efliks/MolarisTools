#!/usr/bin/env python

import sys, os.path
sys.path.append ("/home/mikolaj/devel/MolarisTools")

from MolarisTools import MolarisOutputFile2


molaris = MolarisLogFile2 ("determine_atoms.out")
atoms = (( "MG"  ,  3 ,  "MG"   ,  10. ),
         ( "MG"  ,  4 ,  "MG"   ,  10. ),
         ( "PRX" , -1 ,  "H1'"  ,   3. ),
         ( "PRX" , -1 ,  "H4'"  ,   3. ),
         ( "NUX" , -1 ,  "H5'1" ,   3. ),)
molaris.WriteConstrAtoms (atoms)

angleAtoms = (( "PRX" , -1 ,  "O3'" ),
              ( "NUX" , -1 ,  "PA"  ),
              ( "NUX" , -1 ,  "O3A" ),)
angles = ((angleAtoms, 3., 180.),)
molaris.WriteConstrAngles (angles)
