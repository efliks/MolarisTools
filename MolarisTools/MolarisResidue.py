#-------------------------------------------------------------------------------
# . File      : MolarisResidue.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------

class MolarisResidue (object):
    """A class to extract residue information from Molaris output files."""

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
                chargeTypes = evbAtoms[atomName]
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
if __name__ == "__main__": pass
