#-------------------------------------------------------------------------------
# . File      : MolarisOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine

import exceptions, collections, math


Protein   =  collections.namedtuple ("Protein"  ,  " ebond    ethet     ephi    eitor    evdw     emumu     ehb_pp  ")
Water     =  collections.namedtuple ("Water"    ,  " ebond    ethet                      evdw     emumu     ehb_ww  ")
Prowat    =  collections.namedtuple ("Prowat"   ,  "                                     evdw     emumu     ehb_pw  ")
Long      =  collections.namedtuple ("Long"     ,  " elong                                                          ")
Ac        =  collections.namedtuple ("Ac"       ,  " evd_ac   evd_acw   ehb_ac  ehb_acw  emumuac  emumuacw          ")
Evb       =  collections.namedtuple ("Evb"      ,  " ebond    ethet     ephi             evdw     emumu     eoff  egaschift  eindq   ebulk ")
Induce    =  collections.namedtuple ("Induce"   ,  " eindp    eindw             ")
Const     =  collections.namedtuple ("Const"    ,  " ewatc    eproc     edistc  ")
Langevin  =  collections.namedtuple ("Langevin" ,  " elgvn    evdw_lgv  eborn   ")
Classic   =  collections.namedtuple ("Classic"  ,  " classic  quantum           ")
System    =  collections.namedtuple ("System"   ,  " epot     ekin      etot    ")

MolarisAtom = collections.namedtuple ("molarisAtom", "resSerial resName atomSerial atomName atomType atomCharge x y z connectNames connectSerials")


class MDStep (object):
    """A class to hold energies of an MD step."""

    def __init___ (self):
        for att in ("protein", "water", "prowat", "elong", "ac", "evb", "induce", "const", "langevin", "classic", "system", ):
            setattr (self, att, None)


#  Energies for the system at step          0:
#  ------------------------------------------------------------------------
#  protein - ebond    :      2.57 ethet    :      4.29
#            ephi     :      0.00 eitor    :      0.00
#            evdw     :     -0.24 emumu    :      0.00
#            ehb_pp   :      0.00
#
#  water   - ebond    :    674.31 ethet    :    414.66
#            evdw     :    949.15 emumu    :  -8517.96
#            ehb_ww   :      0.00
#
#  pro-wat - evdw     :     -5.33 emumu    :      0.00
#            ehb_pw   :      0.00
#
#  long    - elong    :     89.62
#
#  ac      - evd_ac   :      0.00 emumuac  :      0.00
#            evd_acw  :      0.00 emumuacw :      0.00
#            ehb_ac   :      0.00
#            ehb_acw  :      0.00
#
#  evb     - ebond    :      0.00 ethet    :      0.00 ephi     :      0.00
#            evdw     :     11.93 emumu    :      0.00 eoff     :      0.00
#            egashift :      0.00 eindq    :      0.00 ebulk    :    -99.57
#
#  induce  - eindp    :      0.00 eindw    :      0.00
#
#  const.  - ewatc    :     27.05 eproc    :      1.45 edistc   :     45.08
#
#  langevin- elgvn    :    -33.50 evdw_lgv :     81.45 eborn    :    -33.07
#
#  classic - epot     :  -6518.24 equantum :   -199.94
#
#  system  - epot     :  -6718.18 ekin     :   2140.90 etot     :  -4577.28
#  _____________________________________________________________________________
class MolarisOutputFile (object):
    """A class for reading output files from Molaris."""

    def __init__ (self, filename="rs_fep.out"):
        """Constructor."""
        self.filename      = filename
        self.currentMDStep = None
        # . fepSteps are the FEP steps (usually 11), each consists of many MD steps (usually 500)
        self.fepSteps      = []
        self._Parse ()


    # . Returns the number of FEP steps (lambda 0 ... 1)
    @property
    def nfepSteps (self):
        return len (self.fepSteps)

    # . Returns the total number of MD steps
    @property
    def nmdSteps (self):
        return sum (map (lambda fepStep: len (fepStep), self.fepSteps))


    def _ReadProtein (self, line, lines):
        toka = line.split ()
        tokb = lines.next ().split ()
        tokc = lines.next ().split ()
        tokd = lines.next ().split ()
        protein = Protein (
                ebond  = float ( toka[4] ) ,
                ethet  = float ( toka[7] ) ,
                ephi   = float ( tokb[2] ) ,
                eitor  = float ( tokb[5] ) ,
                evdw   = float ( tokc[2] ) ,
                emumu  = float ( tokc[5] ) ,
                ehb_pp = float ( tokd[2] ) ,
                          )
        self.currentMDStep.protein = protein


    def _ReadSystem (self, line, lines):
        toka   = line.split ()
        system = System (
                epot = float ( toka[4]  ) ,
                ekin = float ( toka[7]  ) ,
                etot = float ( toka[10] ) ,
                        )
        self.currentMDStep.system = system


    def _ReadClassic (self, line, lines):
        toka     = line.split ()
        energies = Classic (
                classic = float ( toka[4] ) ,
                quantum = float ( toka[7] ) ,
                           )
        self.currentMDStep.classic = energies


    def _Parse (self):
        lines   = iter (open (self.filename).readlines ())
        mdSteps = []
        try:
            while True:
                line = lines.next ()
                if line.startswith (" Energies for the system at step"):
                    self.currentMDStep = MDStep ()
                    while True:
                        line = lines.next ()
                        if   line.startswith ( " protein"  ):
                            self._ReadProtein (line, lines)
                        elif line.startswith ( " water"    ):
                            pass
                        elif line.startswith ( " pro-wat"  ):
                            pass
                        elif line.startswith ( " long"     ):
                            pass
                        elif line.startswith ( " ac"       ):
                            pass
                        elif line.startswith ( " evb"      ):
                            pass
                        elif line.startswith ( " induce"   ):
                            pass
                        elif line.startswith ( " const."   ):
                            pass
                        elif line.startswith ( " langevin" ):
                            pass
                        elif line.startswith ( " classic"  ):
                            self._ReadClassic (line, lines)
                        elif line.startswith ( " system"   ):
                            self._ReadSystem (line, lines)
                            break
                    mdSteps.append (self.currentMDStep)

                elif line.startswith (" Average energies for the system at the step"):
                    self.fepSteps.append (mdSteps)
                    mdSteps = []
        except StopIteration:
            pass


#-------------------------------------------------------------------------------
class MolarisOutputFile2 (object):
    """Alternative class for reading output files from Molaris used in defining constrained atoms."""

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
if __name__ == "__main__": pass
