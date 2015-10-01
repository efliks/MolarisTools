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
                    entry = "%s%10d%10.3f%6s" % (entry, pdbNumber, evbCharge, evbType)
                    pdbNumbers.append (pdbNumber)
                else:
                    entry = "%s%10.3f%6s" % (entry, evbCharge, evbType)
            entry = "%s      # %6.3f%4s    %-4s" % (entry, pdbCharge, pdbType, atomName)
            print ("%sevb_atm%s" % ("    " * tab, entry))

        if verbose:
            for form, charge in zip (forms, totalEVBCharges):
                print ("# Form %d: Total charge of EVB atoms is %.3f" % (form, charge))
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
# . Charges for the RS are now as in Enzymix and calculated with CHELPG
# . Atom types are taken from Han
nameToChargeType = {
 "PA"   :   (( 1.29 ,   "P0")  ,   ( 1.175  ,   "P+")  ,   ( 1.06 ,   "P0"))  ,#
 "O1A"  :   ((-0.85 ,   "O-")  ,   (-0.77   ,   "O-")  ,   (-0.69 ,   "O-"))  ,#
 "O2A"  :   ((-0.85 ,   "O-")  ,   (-0.77   ,   "O-")  ,   (-0.69 ,   "O-"))  ,#
 "O3A"  :   ((-0.64 ,   "Op")  ,   (-0.80   ,   "Op")  ,   (-0.96 ,   "O-"))  ,#
\
 "PB"   :   (( 1.25 ,   "P0")  ,   ( 1.19   ,   "P0")  ,   ( 1.13 ,   "P0"))  ,#
 "O1B"  :   ((-0.82 ,   "O-")  ,   (-0.89   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O2B"  :   ((-0.82 ,   "O-")  ,   (-0.89   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O3B"  :   ((-0.51 ,   "Op")  ,   (-0.505  ,   "Op")  ,   (-0.50 ,   "Op"))  ,#
\
 "PG"   :   (( 1.11 ,   "P0")  ,   ( 1.12   ,   "P0")  ,   ( 1.13 ,   "P0"))  ,#
 "O1G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O2G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
 "O3G"  :   ((-0.88 ,   "O-")  ,   (-0.92   ,   "O-")  ,   (-0.96 ,   "O-"))  ,#
\
 "O5'"  :   ((-0.52 ,   "O0")  ,   (-0.50   ,   "O0")  ,   (-0.48 ,   "O0"))  ,#
\
 "O3'"  :   ((-0.84 ,   "O3")  ,   (-0.705  ,   "Oq")  ,   (-0.57 ,   "Oq"))  ,#
 "C3'"  :   ((-0.13 ,   "C0")  ,   ( 0.115  ,   "C0")  ,   ( 0.36 ,   "C0"))  ,#
 "H3'"  :   (( 0.14 ,   "H0")  ,   ( 0.13   ,   "H0")  ,   ( 0.12 ,   "H0"))  ,#
 "C2'"  :   ((-0.31 ,   "C0")  ,   (-0.26   ,   "C0")  ,   (-0.21 ,   "C0"))  ,#
 "H2'1" :   (( 0.07 ,   "H0")  ,   ( 0.06   ,   "H0")  ,   ( 0.05 ,   "H0"))  ,#
 "H2'2" :   (( 0.07 ,   "H0")  ,   ( 0.06   ,   "H0")  ,   ( 0.05 ,   "H0"))  ,#
}

resonanceForms    = (2, 3) # (1, 2)
verbose           = False
tab               = 2

atomsNucleotide   = ("O5'", "PA", "O1A", "O2A", "O3A", "PB", "O1B", "O2B", "O3B", "PG", "O1G", "O2G", "O3G", )
atomsPrimer       = ("O3'", "C3'", "H3'", "C2'", "H2'1", "H2'2", )
nucleotide        = Residue ("NUX" , 1 , atomsNucleotide , "determine_atoms.out")
primer            = Residue ("PRX" , 2 , atomsPrimer     , "determine_atoms.out")

nucleotide.WriteEVBAtoms (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)
primer.WriteEVBAtoms     (nameToChargeType, resonanceForms, verbose=verbose, tab=tab)

nucleotide.WriteEVBBonds (nameToChargeType, resonanceForms, tab=tab)
primer.WriteEVBBonds     (nameToChargeType, resonanceForms, tab=tab)
