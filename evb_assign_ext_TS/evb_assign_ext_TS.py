#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : evb_assign.py
# . Copyright : USC, Mikolaj J. Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
"""Parse the output file from Molaris and assign atomic charges and types to the EVB atoms."""


class Residue (object):
    def __init__ (self, residueName, residueNumber, listOfAtomNames, logFile):
        lines         =  open (logFile).readlines ()
        self.atoms    =  {}
        foundResidue  =  False
        for line in lines:
            if not foundResidue:
                if line.count ("atom list for residue"):
                    tokens    = line.split ()
                    residue   = tokens[4]
                    resNumber, resName = residue.split ("_")
                    resNumber = int (resNumber)
                    resName   = resName.replace (",", "")
                    if resNumber == residueNumber and resName == residueName:
                        foundResidue = True
            else:
                if line.count ("Total charge"):
                    break
                tokens = line.split ()
                if len (tokens) > 0:
                    if tokens[0].isdigit ():
                        (atomNumber, atomName, atomType), atomCharge = tokens[:3], tokens[6]
                        atomNumber     = int (atomNumber)
                        atomCharge     = float (atomCharge)
                        atoms          = tokens[7:]
                        connectNames   = []
                        connectNumbers = []
                        for atom in atoms:
                            if atom.isdigit ():
                                connectNumbers.append (int (atom))
                            else:
                                connectNames.append (atom)
                        if atomName in listOfAtomNames:
                            self.atoms[atomName] = [atomCharge, atomType, atomNumber, connectNames, connectNumbers]
    

    def WriteEVBAtoms (self, evbAtoms, forms, verbose=False, tab=0):
        totalEVBCharges = [0.] * len (forms)
        pdbNumbers      = []
        for atomName, atomData in self.atoms.iteritems ():
            if not evbAtoms.has_key (atomName):
                print ("# EVB atom %s not present" % atomName)
                continue
            chargeTypes = evbAtoms[atomName]
            pdbCharge, pdbType, pdbNumber, connectNames, connectNumbers = atomData
            entry = ""
            for fi, form in enumerate (forms, 0):
                evbCharge, evbType   = chargeTypes[form - 1]
                totalEVBCharges[fi] += evbCharge
                if fi < 1:
                    entry = "%s%10d%10.4f%6s" % (entry, pdbNumber, evbCharge, evbType)
                    pdbNumbers.append (pdbNumber)
                else:
                    entry = "%s%10.4f%6s" % (entry, evbCharge, evbType)
            entry = "%s      # %10.4f%4s    %-4s" % (entry, pdbCharge, pdbType, atomName)
            print ("%sevb_atm%s" % ("    " * tab, entry))

        if verbose:
            for form, charge in zip (forms, totalEVBCharges):
                print ("# Form %d: Total charge of EVB atoms is %.4f" % (form, charge))
            print ("# %s" % (" ".join (map (str, pdbNumbers))))


    def WriteEVBBonds (self, evbAtoms, forms, tab=0):
        bonds       = []
        typeChanges = {}

        for atomName, atomData in self.atoms.iteritems ():
            isEVBAtom = evbAtoms.has_key (atomName)
            if isEVBAtom:
                chargeTypes = nameToChargeType[atomName]
                pdbCharge, pdbType, atomNumber, connectNames, connectNumbers = atomData

                # . Determine if the atom changes its type between resonance forms
                typeChanges[atomName] = False
                prevType        = chargeTypes[forms[0] - 1][1]
                for form in forms[1:]:
                    evbCharge, evbType = chargeTypes[form - 1]
                    if evbType != prevType:
                        typeChanges[atomName] = True
                        break
                    prevType = evbType

                for otherAtomNumber, otherAtomName in zip (connectNumbers, connectNames):
                    # . The other atom has to be an EVB atom as well
                    if evbAtoms.has_key (otherAtomName):

                        pair = [(atomNumber, atomName), (otherAtomNumber, otherAtomName)]
                        if pair not in bonds:
                            pair.reverse ()
                            if pair not in bonds:
                                bonds.append (pair)
        for pair in bonds:
            (aNumber, aName), (bNumber, bName) = pair
            if typeChanges.has_key (aName):
                aChanges = typeChanges[aName]
            else:
                aChanges = False
            if typeChanges.has_key (bName):
                bChanges = typeChanges[bName]
            else:
                bChanges = False
            print ("%sevb_bnd   0    %6d%6d   # %-6s  %-6s%s" % ("    " * tab, aNumber, bNumber, aName, bName, "  (*)" if any ((aChanges, bChanges)) else ""))


#===============================================================================
# . Main program
#===============================================================================
# . Charges for the RS are now as in Enzymix and calculated with CHELPG.
# . Atom types are taken from Han.
# . TS1, INT and TS2 share the same atom types.
#                  RS                   TS1                 INT                 TS2                 PS
nameToChargeType = {
"PA"   :   (( 1.24  , "P0")  ,  ( 1.195  , "P+") ,  ( 1.15   , "P+") ,   ( 1.105  , "P+") ,   ( 1.06  , "P0")) ,
"O1A"  :   ((-0.71  , "O-")  ,  (-0.705  , "O-") ,  (-0.7    , "O-") ,   (-0.695  , "O-") ,   (-0.69  , "O-")) ,
"O2A"  :   ((-0.71  , "O-")  ,  (-0.705  , "O-") ,  (-0.7    , "O-") ,   (-0.695  , "O-") ,   (-0.69  , "O-")) ,
"O3A"  :   ((-0.41  , "Op")  ,  (-0.48   , "Op") ,  (-0.55   , "Op") ,   (-0.62   , "Op") ,   (-0.69  , "O-")) ,
"O5'"  :   ((-0.41  , "O0")  ,  (-0.4275 , "O0") ,  (-0.445  , "O0") ,   (-0.4625 , "O0") ,   (-0.48  , "O0")) ,
\
"O3'"  :   ((-0.90  , "O3")  ,  (-0.7875 , "Oq") ,  (-0.675  , "Oq") ,   (-0.5625 , "Oq") ,   (-0.45  , "Oq")) ,
"C3'"  :   ((-0.25  , "C0")  ,  (-0.225  , "C0") ,  (-0.2    , "C0") ,   (-0.175  , "C0") ,   (-0.15  , "C0")) ,
"H3'1" :   (( 0.05  , "H0")  ,  ( 0.045  , "H0") ,  ( 0.04   , "H0") ,   ( 0.035  , "H0") ,   ( 0.03  , "H0")) ,
"H3'2" :   (( 0.05  , "H0")  ,  ( 0.045  , "H0") ,  ( 0.04   , "H0") ,   ( 0.035  , "H0") ,   ( 0.03  , "H0")) ,
"H3'3" :   (( 0.05  , "H0")  ,  ( 0.045  , "H0") ,  ( 0.04   , "H0") ,   ( 0.035  , "H0") ,   ( 0.03  , "H0")) ,
}

stateToNumber = { "RS"  : 1 ,
                  "TS1" : 2 ,
                  "INT" : 3 ,
                  "TS2" : 4 ,
                  "PS"  : 5 , }

resonanceForms    = (stateToNumber["RS"], stateToNumber["TS1"])
verbose           = False
tab               = 2

atomsNucleotide   = ("PA", "O1A", "O2A", "O3A", "O5'", )
atomsPrimer       = ("O3'", "C3'" , "H3'1", "H3'2", "H3'3", )
nucleotide        = Residue ("NUS" , 1 , atomsNucleotide , "../determine_atoms.out")
primer            = Residue ("PRS" , 2 , atomsPrimer     , "../determine_atoms.out")

nucleotide.WriteEVBAtoms (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)
primer.WriteEVBAtoms     (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)

nucleotide.WriteEVBBonds (nameToChargeType, resonanceForms, tab=tab)
primer.WriteEVBBonds     (nameToChargeType, resonanceForms, tab=tab)
