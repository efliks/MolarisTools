#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : molaris_constr.py
# . Copyright : USC, Mikolaj J. Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------

import collections
import exceptions

MolarisAtom = collections.namedtuple ("molarisAtom", "resSerial resName atomSerial atomName atomType atomCharge x y z connectNames connectSerials")


class MolarisLogFile (object):
    def __init__ (self, logFile):
        lines         =  open (logFile).readlines ()
        self.atoms    =  []
        foundResidue  =  False
        for line in lines:
            if not foundResidue:
                if line.count ("atom list for residue"):
                    tokens       = line.split ()
                    residue      = tokens[4]
                    resSerial, resName = residue.split ("_")
                    resSerial    = int (resSerial)
                    resName      = resName.replace (",", "")
                    foundResidue = True
            else:
                if line.count ("Total charge"):
                    foundResidue = False
                else:
                    tokens = line.split ()
                    if len (tokens) > 0:
                        if tokens[0].isdigit ():
                            (atomSerial, atomName, atomType), (x, y, z), atomCharge = tokens[:3], map (float, tokens[3:6]), float (tokens[6])
                            atomSerial     = int (atomSerial)
                            atoms          = tokens[7:]
                            connectNames   = []
                            connectSerials = []
                            for atom in atoms:
                                if atom.isdigit ():
                                    connectSerials.append (int (atom))
                                else:
                                    connectNames.append (atom)
                            molaris = MolarisAtom (x=x, y=y, z=z, resName=resName, resSerial=resSerial, atomName=atomName, atomType=atomType, atomSerial=atomSerial, atomCharge=atomCharge, connectNames=connectNames, connectSerials=connectSerials,)
                            self.atoms.append (molaris)


    def _FindAtom (self, resName, resSerial, atomName):
        found = False
        if resSerial > 0:
            for atom in self.atoms:
                if atom.resSerial == resSerial and atom.resName == resName and atom.atomName == atomName:
                    found = True
                    break
        else:
            # . If the serial number is lower than zero, the residue has to have a unique name
            for atom in self.atoms:
                if atom.resName == resName and atom.atomName == atomName:
                    found = True
                    break
        if not found:
            raise exceptions.StandardError ("Cannot find atom %s %d %s" % (resName, resSerial, atomName))
        return atom


    def WriteConstrAtoms (self, atoms):
        for resName, resSerial, atomName, force in atoms:
            atom = self._FindAtom (resName, resSerial, atomName)
            print ("constraint_post  %4d     %4.1f   %4.1f   %4.1f    %8.3f  %8.3f  %8.3f   # %s" % (atom.atomSerial, force, force, force, atom.x, atom.y, atom.z, atom.atomName))


    def WriteConstrAngles (self, angles):
        for atoms, force, equilAngle in angles:
            found = []
            for resName, resSerial, atomName in atoms:
                atom = self._FindAtom (resName, resSerial, atomName)
                found.append (atom)
            print ("constraint_ang  %4d   %4d   %4d   %4.1f   %8.1f   # %4s %4s %4s" % (found[0].atomSerial, found[1].atomSerial, found[2].atomSerial, force, equilAngle, found[0].atomName, found[1].atomName, found[2].atomName))


#===============================================================================
# . Main program
#===============================================================================
molaris = MolarisLogFile ("determine_atoms.out")

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
