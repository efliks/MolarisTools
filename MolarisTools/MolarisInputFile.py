#-------------------------------------------------------------------------------
# . File      : MolarisInputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine
import  collections, exceptions


EVBAtom  = collections.namedtuple ("EVBAtom"  , "atomSerial  atomTypeI  atomChargeI  atomTypeII  atomChargeII")
EVBBond  = collections.namedtuple ("EVBBond"  , "states  bonda  bondb")



class MolarisInputFile (object):
    """A class to represent a Molaris input file."""

    def __init__ (self, filename):
        """Constructor."""
        self.filename = filename
        self._Parse ()


    def _GetLineWithComment (self, data, beginComment="#"):
        line      = data.next ()
        position  = line.find (beginComment)
        if position > -1:
            text      = line[             : position].strip ()
            comment   = line[position + 1 :         ].strip ()
        else:
            text      = line
            comment   = ""
        return (text, comment)


    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0

    @property
    def nbonds (self):
        if hasattr (self, "bonds"):
            return len (self.bonds)
        return 0

    @property
    def charge (self):
        if hasattr (self, "atoms"):
            charge = 0.
            for atom in self.atoms:
                charge += atom.atomChargeI
        return 0.

    @property
    def chargeII (self):
        if hasattr (self, "atoms"):
            charge = 0.
            for atom in self.atoms:
                charge += atom.atomChargeII
        return 0.

    @property
    def types (self):
        if hasattr (self, "atoms"):
            types = []
            for atom in self.atoms:
                if atom.atomTypeI not in types:
                    types.append (atom.atomTypeI)
            types.sort ()
            return types
        return []

    @property
    def typesII (self):
        if hasattr (self, "atoms"):
            types = []
            for atom in self.atoms:
                if atom.atomTypeII not in types:
                    types.append (atom.atomTypeII)
            types.sort ()
            return types
        return []

    @property
    def ntypes (self):
        return len (self.types)

    @property
    def ntypesII (self):
        return len (self.typesII)


    def _Parse (self):
        lines      = open (self.filename)
        self.atoms = []
        self.bonds = []
        try:
            while True:
                line, comment = self._GetLineWithComment (lines)

                # . Read EVB states, assume there are only two of them
                # evb_atm         6   -0.3000    C0   -0.3000    C0
                # evb_atm         5   -0.4500    O0   -0.4500    O0
                # (...)
                if line.startswith ("evb_atm"):
                    foo, atomSerial, atomChargeI, atomTypeI, atomChargeII, atomTypeII = TokenizeLine (line, converters=[None, int, float, None, float, None])
                    atom = EVBAtom (
                        atomSerial   = atomSerial    , 
                        atomChargeI  = atomChargeI   , 
                        atomTypeI    = atomTypeI     , 
                        atomChargeII = atomChargeII  , 
                        atomTypeII   = atomTypeII    ,
                        )
                    self.atoms.append (atom)

                # . Read EVB bonds
                #        evb_bnd   0         5     6   # O5'     C5'   
                #        evb_bnd   0         7     6   # H5'1    C5'   
                # (...)
                elif line.startswith ("evb_bnd"):
                    foo, states, bonda, bondb = TokenizeLine (line, converters=[None, int, int, int])
                    bond = EVBBond (
                        states  =   states  ,
                        bonda   =   bonda   ,
                        bondb   =   bondb   ,
                        )
                    self.bonds.append (bond)


#                if line.count ("constraint_pair"):
#                    foo, aserial, bserial, forceConst, equilDist = TokenizeLine (line, converters=[None, int, int, float, float])
#                    pair  = (aserial, bserial)
#                    found = False
#                    if pair == self.pair:
#                        found = True
#                    else:
#                        pair = (bserial, aserial)
#                        if pair == self.pair:
#                            found = True
#                    if found:
#                        equil = [forceConst, equilDist]
#                        break
        except StopIteration:
            pass
        # . Close the file
        lines.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
