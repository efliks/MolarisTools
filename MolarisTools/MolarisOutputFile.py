#-------------------------------------------------------------------------------
# . File      : MolarisOutputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
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
Evb       =  collections.namedtuple ("Evb"      ,  " ebond    ethet     ephi             evdw     emumu     eoff  egashift   eindq   ebulk ")
Induce    =  collections.namedtuple ("Induce"   ,  " eindp    eindw             ")
Const     =  collections.namedtuple ("Const"    ,  " ewatc    eproc     edistc  ")
Langevin  =  collections.namedtuple ("Langevin" ,  " elgvn    evdw_lgv  eborn   ")
Classic   =  collections.namedtuple ("Classic"  ,  " classic  quantum           ")
System    =  collections.namedtuple ("System"   ,  " epot     ekin      etot    ")

MolarisAtom = collections.namedtuple ("molarisAtom", "resSerial resName atomSerial atomName atomType atomCharge x y z connectNames connectSerials")


Atom      = collections.namedtuple ("Atom"    ,   "label  serial  charge  bonds  x  y  z")
Residue   = collections.namedtuple ("Residue" ,   "label  serial  atoms")



class MDStep (object):
    """A class to hold energies of an MD step."""

    def __init___ (self):
        for att in ("protein", "water", "prowat", "elong", "ac", "evb", "induce", "const", "langevin", "classic", "system", ):
            setattr (self, att, None)


class MolarisOutputFile (object):
    """A class for reading output files from Molaris."""

    def __init__ (self, filename="rs_fep.out", logging=False):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    # . Returns the number of FEP steps (lambda 0 ... 1)
    @property
    def nfepSteps (self):
        if hasattr (self, "fepSteps"):
            return len (self.fepSteps)
        return 0

    # . Returns the total number of MD steps
    @property
    def nmdSteps (self):
        if hasattr (self, "fepSteps"):
            nmdSteps = 0
            for fepStep in self.fepSteps:
               nmdSteps += len (fepStep)
            return nmdSteps
        return 0

    # . Returns the number of residues (useful for determine_atoms type of script)
    @property
    def nresidues (self):
        if hasattr (self, "residues"):
            return len (self.residues)
        return 0

    # . Returns the total number of atoms from all residues
    @property
    def natoms (self):
        if hasattr (self, "residues"):
            total = 0
            for residue in self.residues:
                total += len (residue.atoms)
            return total
        return 0


    def _Parse (self, logging=False):
        # . fepSteps are the FEP steps (usually 11), each consisting of many MD steps (usually 500)
        fepSteps      = []
        mdSteps       = []
        currentMDStep = None
        residues      = []
        lines         = open (self.filename)
        try:
            while True:
                line = lines.next ()

                #  Energies for the system at step          0:
                #  ------------------------------------------------------------------------
                if line.startswith (" Energies for the system at step"):
                    currentMDStep = MDStep ()
                    while True:
                        line = lines.next ()

                        #  protein - ebond    :      2.57 ethet    :      4.29
                        #            ephi     :      0.00 eitor    :      0.00
                        #            evdw     :     -0.24 emumu    :      0.00
                        #            ehb_pp   :      0.00
                        #
                        if   line.startswith ( " protein"  ):
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
                                    ehb_pp = float ( tokd[2] ) ,)
                            currentMDStep.protein = protein


                        #  water   - ebond    :    674.31 ethet    :    414.66
                        #            evdw     :    949.15 emumu    :  -8517.96
                        #            ehb_ww   :      0.00
                        #
                        elif line.startswith ( " water"    ):
                            pass
                        #  pro-wat - evdw     :     -5.33 emumu    :      0.00
                        #            ehb_pw   :      0.00
                        #
                        elif line.startswith ( " pro-wat"  ):
                            pass
                        #  long    - elong    :     89.62
                        #
                        elif line.startswith ( " long"     ):
                            pass
                        #  ac      - evd_ac   :      0.00 emumuac  :      0.00
                        #            evd_acw  :      0.00 emumuacw :      0.00
                        #            ehb_ac   :      0.00
                        #            ehb_acw  :      0.00
                        #
                        elif line.startswith ( " ac"       ):
                            pass
                        #  evb     - ebond    :      0.00 ethet    :      0.00 ephi     :      0.00
                        #            evdw     :     11.93 emumu    :      0.00 eoff     :      0.00
                        #            egashift :      0.00 eindq    :      0.00 ebulk    :    -99.57
                        #
                        elif line.startswith ( " evb"      ):
                            toka = line.split ()
                            tokb = lines.next ().split ()
                            tokc = lines.next ().split ()
                            evb  = Evb (
                                ebond    = float (toka[4] )    ,
                                ethet    = float (toka[7] )    ,
                                ephi     = float (toka[10])    ,
                                evdw     = float (tokb[2] )    ,
                                emumu    = float (tokb[5] )    ,
                                eoff     = float (tokb[8] )    ,
                                egashift = float (tokc[2] )    ,
                                eindq    = float (tokc[5] )    ,
                                ebulk    = float (tokc[8] )    ,)
                            currentMDStep.evb = evb

                        #  induce  - eindp    :      0.00 eindw    :      0.00
                        #
                        elif line.startswith ( " induce"   ):
                            pass
                        #  const.  - ewatc    :     27.05 eproc    :      1.45 edistc   :     45.08
                        #
                        elif line.startswith ( " const."   ):
                            pass
                        #  langevin- elgvn    :    -33.50 evdw_lgv :     81.45 eborn    :    -33.07
                        #
                        elif line.startswith ( " langevin" ):
                            pass
                        #  classic - epot     :  -6518.24 equantum :   -199.94
                        #
                        elif line.startswith ( " classic"  ):
                            toka     = line.split ()
                            energies = Classic (
                                    classic = float ( toka[4] ) ,
                                    quantum = float ( toka[7] ) ,)
                            currentMDStep.classic = energies


                        #  system  - epot     :  -6718.18 ekin     :   2140.90 etot     :  -4577.28
                        #  _____________________________________________________________________________
                        elif line.startswith ( " system"   ):
                            toka   = line.split ()
                            system = System (
                                    epot = float ( toka[4]  ) ,
                                    ekin = float ( toka[7]  ) ,
                                    etot = float ( toka[10] ) ,)
                            currentMDStep.system = system
                            break
                    mdSteps.append (currentMDStep)
                    if logging:
                        nsteps = len (mdSteps)
                        print ("Read energies for %d MD step%s." % (nsteps, "" if nsteps < 2 else "s"))

                elif line.startswith (" Average energies for the system at the step"):
                    fepSteps.append (mdSteps)
                    mdSteps = []
                    if logging:
                        nfep = len (fepSteps)
                        print ("Read FEP step %d." % nfep)



                # atom list for residue:     2_WAT,    # of atoms in this residue:   3
                #
                # number  name  type      x          y          z      charge      atoms bonded(name)      atoms bonded(number)
                # ------  ----  ----   -------    -------    -------   ------   ------------------------ ------------------------
                elif line.count ("atom list for residue"):
                    tokens    = line.split ()
                    residue   = tokens[4]
                    resSerial, resLabel = residue.split ("_")
                    resSerial = int (resSerial)
                    resLabel  = resLabel.replace (",", "")
                    # . Skip a few lines
                    for i in range (3):
                        lines.next ()
                    # . Read atoms
                    #    2    OH     O2     -0.087     -0.022      2.081   -0.800   H1   H2                      3     4
                    #    3    H1     H2     -0.139     -0.807      2.652    0.400   OH                           2
                    #    4    H2     H2     -0.056      0.751      2.671    0.400   OH                           2
                    #
                    # Total charge of this residue:     0.000
                    atoms     = []
                    while True:
                        line   = lines.next ()
                        if line.count ("Total charge"):
                            break
                        tokens = line.split ()
                        if len (tokens) > 0:
                            if tokens[0].isdigit ():
                                (atomSerial, atomLabel, atomType), atomCharge = tokens[:3], tokens[6]
                                atomSerial     = int (atomSerial)
                                atomCharge     = float (atomCharge)
                                x, y, z        = map (float, tokens[3:6])
                                # . Read atoms the current atom is connected to
                                bondAtoms      = tokens[7:]
                                bondLabels     = []
                                bondSerials    = []
                                for atom in bondAtoms:
                                    if atom.isdigit ():
                                        bondSerials.append (int (atom))
                                    else:
                                        bondLabels.append (atom)
                                # . Prepare bonded atoms
                                bonds = []
                                for serial, label in zip (bondSerials, bondLabels):
                                    bondedAtom = (serial, label)
                                    bonds.append (bondedAtom)
                                # . Add a new atom
                                atom = Atom (
                                    label   =   atomLabel   ,
                                    serial  =   atomSerial  ,
                                    charge  =   atomCharge  ,
                                    bonds   =   bonds       ,
                                    x       =   x           ,
                                    y       =   y           ,
                                    z       =   z           ,
                                    )
                                atoms.append (atom)
                    # . Add a new residue
                    residue = Residue (
                        serial  =   resSerial   ,
                        label   =   resLabel    ,
                        atoms   =   atoms       ,
                        )
                    residues.append (residue)
                    if logging:
                        natoms = len (atoms)
                        print ("Found residue %s-%d with %d atom%s." % (resLabel, resSerial, natoms, "" if natoms < 2 else "s"))


                # Classical forces which are not calculated in qm:
                #   evb_atom     fx        fy        fz
                #         6     2.792   -10.205     9.635
                #         5    -7.012    14.479   -21.646
                # (...)
                elif line.count ("Classical forces which are not calculated in qm"):
                    for i in range (2):
                        line = lines.next ()
                    forces = []
                    while line != "\n":
                        tokens = TokenizeLine (line, converters=[int, float, float, float])
                        evbSerial, fx, fy, fz = tokens
                        force  = (fx, fy, fz)
                        forces.append (force)
                        line   = lines.next ()
                    #if logging:
                    #    nforces = len (forces)
                    #    print ("Read classical forces for %d atoms." % nforces)


                #  Forces(classical+qm) and Charges will be used for dynamics:
                #  EVB_atom	 x         y         z        fx        fy       fz      crg
                #         6     3.666     6.485    12.973    -8.915    -1.595    -4.922   0.272
                #         5     4.872     6.423    13.671     1.443     1.123     0.273  -0.750
                # (...)
                elif line.count ("Forces(classical+qm) and Charges will be used for dynamics:"):
                    for i in range (2):
                        line = lines.next ()
                    forcesQMMM = []
                    while line != "\n":
                        tokens = TokenizeLine (line, converters=[int, float, float, float, float, float, float, float])
                        evbSerial, x, y, z, fx, fy, fz, charge = tokens
                        force  = (fx, fy, fz)
                        forcesQMMM.append (force)
                        line   = lines.next ()
                    #if logging:
                    #    nforces = len (forces)
                    #    print ("Read QM/MM forces for %d atoms." % nforces)
                    if logging:
                        for forceClassic, forceQMMM in zip (forces, forcesQMMM):
                            fx, fy, fz = forceClassic
                            qmmmFx, qmmmFy, qmmmFz = forceQMMM
                            dfx = qmmmFx - fx
                            dfy = qmmmFy - fy
                            dfz = qmmmFz - fz
                            print ("%8.3f  %8.3f  %8.3f      %8.3f  %8.3f  %8.3f      %8.3f  %8.3f  %8.3f" % (fx, fy, fz, qmmmFx, qmmmFy, qmmmFz, dfx, dfy, dfz))
                        print ("======================")

        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . Finish up
        if fepSteps != []:
            self.fepSteps = fepSteps
        if currentMDStep != None:
            self.currentMDStep = currentMDStep
        if residues != []:
            self.residues = residues


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
