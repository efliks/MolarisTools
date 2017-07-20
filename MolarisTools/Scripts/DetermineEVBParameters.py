#-------------------------------------------------------------------------------
# . File      : DetermineEVBParameters.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, math

from MolarisTools.Units    import DEFAULT_EVB_LIB
from MolarisTools.Parser   import MolarisInputFile, EVBDatFile
from MolarisTools.Library  import EVBLibrary, EVBMorseAtom, EVBMorsePair


def DetermineEVBParameters (filenameInput="heat_template.inp", filenameDat=os.path.join ("evb_heat_01", "evb.dat"), filenameEVBLibrary=DEFAULT_EVB_LIB, state=1, logging=True):
    """Get EVB parameters for a system.

    For now, only bonding parameters are returned."""
    # . Dat file contains lists of parameters used by Molaris for a particular system
    dat      = EVBDatFile       (filenameDat        , logging=logging)
    # . Input file is used to relate atom serial numbers to their types
    mif      = MolarisInputFile (filenameInput      , logging=logging)
    # . Library is used to pick up the currently used parameters
    library  = EVBLibrary       (filenameEVBLibrary , logging=logging)

    convert  = {}
    for (i, atom) in enumerate (mif.states[(state - 1)], 1):
        # convert[atom.serial] = atom.atype
        # . Molaris does not seem to use the serial numbers of atoms, 
        #       but renumbers them in the order as they appear in the input file
        convert[i] = atom.atype
    # . Generate a unique list of atom types in bonds
    bonds    = []
    for bond in dat.bonds:
        if bond.exist[(state - 1)]:
            (seriala , serialb) = bond.serials
            try:
                (typea   , typeb  ) = (convert[seriala], convert[serialb])
                pair = (typea, typeb)
                riap = (typeb, typea)
                if not ((pair in bonds) or (riap in bonds)):
                    bonds.append (pair)
            except:
                if logging:
                    print ("# . Warning: Bond (%d, %d, %d) involves non-EVB atoms" % (seriala, serialb))
    bonds.sort ()
    # . Generate a unique list of atom types in angles
    angles   = []
    for angle in dat.angles:
        if angle.exist[(state - 1)]:
            (seriala , serialb , serialc) = angle.serials
            try:
                (typea   , typeb   , typec  ) = (convert[seriala], convert[serialb], convert[serialc])
                triplet = (typea, typeb, typec)
                telpirt = (typec, typeb, typea)
                # if not ((triplet in angles) or (telpirt in angles)):
                found   = False
                for (storeda, storedb, storedc) in angles:
                    if storedb == typeb:
                        found = True
                        break
                if not found:
                    angles.append (triplet)
            except:
                if logging:
                    print ("# . Warning: Angle (%d, %d, %d) involves non-EVB atoms" % (seriala, serialb, serialc))
    angles.sort ()
    # . Generate a unique list of atom types in dihedral angles
    torsions = []
    for torsion in dat.torsions:
        if torsion.exist[(state - 1)]:
            (seriala , serialb , serialc , seriald) = torsion.serials
            try:
                (typea   , typeb   , typec   , typed  ) = (convert[seriala], convert[serialb], convert[serialc], convert[seriald])
                quartet = (typea, typeb, typec, typed)
                tetrauq = (typed, typec, typeb, typea)
                # if not ((quartet in torsions) or (tetrauq in torsions)):
                found   = False
                for (storeda, storedb, storedc, storedd) in torsions:
                    if ((storedb == typeb) and (storedc == typec)) or ((storedb == typec) and (storedc == typeb)):
                        found = True
                        break
                if not found:
                    torsions.append (quartet)
            except:
                if logging:
                    print ("# . Warning: Torsion angle (%d, %d, %d, %d) involves non-EVB atoms" % (seriala, serialb, serialc, seriald))
    torsions.sort ()

    # . Collect bond parameters from the EVB library
    parBonds = []
    for bond in bonds:
        (typea, typeb) = bond
        parameter = library.GetBond (typea, typeb)
        if parameter:
            if not isinstance (parameter, EVBMorsePair):
                # . Convert individual EVB atoms to an EVB pair
                if logging:
                    print ("# . Applying combination rules for atom types (%s, %s)" % (typea, typeb))
                (morsea, morseb) = parameter
                morseD  = math.sqrt (morsea.morseD * morseb.morseD)
                r0      = morsea.radius + morseb.radius
                pair    = EVBMorsePair (
                    typea           =  typea          ,
                    typeb           =  typeb          ,
                    morseAB         =  morseD         ,
                    rab             =  r0             ,
                    beta            =  _DEFAULT_BETA  ,
                    forceHarmonic   =  400.           ,
                    radiusHarmonic  =    1.4          , )
                parameter = pair
            parBonds.append (parameter)
        else:
            if logging:
                print ("# . Warning: Bond parameters for atom types (%s, %s) not found" % (typea, typeb))
    # . Collect angle parameters from the EVB library
    parAngles = []
    for angle in angles:
        (typea, typeb, typec) = angle
        parameter = library.GetAngle (typea, typeb, typec)
        if parameter:
            parAngles.append (parameter)
        else:
            if logging:
                print ("# . Warning: Angle parameters for atom types (%s, %s, %s) not found" % (typea, typeb, typec))
    # . Collect dihedral parameters from the EVB library
    parTorsions = []
    for torsion in torsions:
        (typea, typeb, typec, typed) = torsion
        parameter = library.GetTorsion (typeb, typec)
        if parameter:
            parTorsions.append (parameter)
        else:
            if logging:
                print ("# . Warning: Torsion parameters for atom types (%s, %s, %s, %s) not found" % (typea, typeb, typec, typed))

    if logging:
        print ("\n-- EVB parameters for bonds --")
        for bond    in parBonds:
            #  EVBMorsePair  = collections.namedtuple ("EVBMorsePair"  ,  "typea  typeb  morseAB  rab  beta  forceHarmonic  radiusHarmonic")
            print ("%2s    %2s    %8.2f    %8.2f    %8.1f" % (bond.typea, bond.typeb, bond.morseAB, bond.rab, bond.beta))
        print ("\n-- EVB parameters for angles --")
        for angle   in parAngles:
            #  EVBAngle      = collections.namedtuple ("EVBAngle"      ,  "evbType  force  angle0  foo  bar")
            print ("%2s    %8.2f    %8.2f" % (angle.evbType, angle.force, angle.angle0))
        print ("\n-- EVB parameters for torsions --")
        for torsion in parTorsions:
            #  EVBTorsion    = collections.namedtuple ("EVBTorsion"    ,  "typea  typeb  force  periodicity  phase")
            print ("%2s    %2s    %8.2f    %8.2f    %8.2f" % (torsion.typea, torsion.typeb, torsion.force, torsion.periodicity, torsion.phase))
    return (parBonds, parAngles, parTorsions)


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
