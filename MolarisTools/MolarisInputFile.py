#-------------------------------------------------------------------------------
# . File      : MolarisInputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Units     import *
from   Utilities import TokenizeLine
import  collections, exceptions


EVBAtom  = collections.namedtuple ("EVBAtom"  , "serial  atype  charge  comment")
EVBBond  = collections.namedtuple ("EVBBond"  , "states  seriala  serialb  comment")

_MODULE_LABEL = "MolarisInput"


class MolarisInputFile (object):
    """A class to represent a Molaris input file."""

    def __init__ (self, filename, logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    def GetPairs (self, state=1):
        """Get a pair (label, type) for each EVB atom.

        Atoms must have comments generated by GenerateEVBList, since these comments contain atom labels:

        evb_atm     1     1.15    P0        1.15    P0    #   0.700    P4    PG    A
        evb_atm     2    -0.76    O-       -0.76    O-    #  -0.900    O3   O1G    A
        (...)"""
        pairs = []
        atoms = self.states[state - 1]
        for atom in atoms:
            (enzymixCharge, enzymixType, label, groupLabel) = TokenizeLine (atom.comment, converters=[float, None, None, None])
            pair = (label, atom.atype)
            pairs.append (pair)
        return pairs


    def WriteStates (self, selection=0, shift=2, comments=True):
        """Write a table of statesi.

        If selection=1 or selection=2, write state I or state II twice, respectively."""
        if hasattr (self, "states"):
            assign = {
                0   :   (self.states[0], self.states[1])   ,
                1   :   (self.states[0], self.states[0])   ,
                2   :   (self.states[1], self.states[1])   , }
            stateI, stateII = assign[selection]
            for atomI, atomII in zip (stateI, stateII):
                comment = ""
                if comments:
                    if atomI.comment:
                        comment = "    #%s" % atomI.comment
                print ("%sevb_atm   %4d     %5.2f   %2s          %5.2f   %2s%s" % ("    " * shift, atomI.serial, atomI.charge, atomI.atype, atomII.charge, atomII.atype, comment))


    def WriteBonds (self, shift=2, comments=True):
        """Write a table of bonds."""
        if hasattr (self, "bonds"):
            for bond in self.bonds:
                comment = ""
                if comments:
                    if bond.comment:
                        comment = "    #%s" % bond.comment
                print ("%sevb_bnd%4d  %4d  %4d%s" % ("    " * shift, bond.states, bond.seriala, bond.serialb, comment))


    def _GetLineWithComment (self, data, beginComment="#"):
        line      = data.next ()
        position  = line.find (beginComment)
        if position > -1:
            text      = line[             : position].strip ()
            comment   = line[position + 1 :         ][:-1]
        else:
            text      = line
            comment   = ""
        return (text, comment)


    @property
    def natoms (self):
        if hasattr (self, "states"):
            return len (self.states[0])
        return 0

    @property
    def nbonds (self):
        if hasattr (self, "bonds"):
            return len (self.bonds)
        return 0

    @property
    def charges (self):
        if hasattr (self, "states"):
            charge, chargeOther = (0., 0.)
            for (atom, atomOther) in zip (self.states[0], self.states[1]):
                charge += atom.charge
                chargeOther += atomOther.charge
            return (charge, chargeOther)
        return (0., 0.)

    @property
    def types (self):
        if hasattr (self, "states"):
            types, typesOther = ([], [])
            for (atom, atomOther) in zip (self.states[0], self.states[1]):
                if atom.atype not in types:
                    types.append (atom.atype)
                if atomOther.atype not in typesOther:
                    typesOther.append (atomOther.atype)
            types.sort ()
            typesOther.sort ()
            return (types, typesOther)
        return ([], [])

    @property
    def ntypes (self):
        (ntypes, ntypesOther) = map (len, self.types)
        return (ntypes, ntypesOther)

    @property
    def nstates (self):
        return 2


    def _Parse (self, logging):
        lines       = open (self.filename)
        stateI      = []
        stateII     = []
        self.bonds  = []
        try:
            while True:
                line, comment = self._GetLineWithComment (lines)

                # . Read EVB states, assume there are only two of them
                # evb_atm         6   -0.3000    C0   -0.3000    C0
                # evb_atm         5   -0.4500    O0   -0.4500    O0
                # (...)
                if line.startswith ("evb_atm"):
                    foo, serial, chargeI, atypeI, chargeII, atypeII = TokenizeLine (line, converters=[None, int, float, None, float, None])
                    atomState = EVBAtom (
                        atype   =   atypeI   ,
                        serial  =   serial   ,
                        charge  =   chargeI  ,
                        comment =   comment  , )
                    stateI.append (atomState)
                    
                    atomState = EVBAtom (
                        atype   =   atypeII   ,
                        serial  =   serial    ,
                        charge  =   chargeII  ,
                        comment =   comment   , )
                    stateII.append (atomState)


                # . Read EVB bonds
                #        evb_bnd   0         5     6   # O5'     C5'   
                #        evb_bnd   0         7     6   # H5'1    C5'   
                # (...)
                elif line.startswith ("evb_bnd"):
                    foo, states, seriala, serialb = TokenizeLine (line, converters=[None, int, int, int])
                    bond = EVBBond (
                        states    =   states   ,
                        seriala   =   seriala  ,
                        serialb   =   serialb  ,
                        comment   =   comment  ,
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
        # . Finalize
        self.states = [stateI, stateII]
        if logging:
            print ("# . %s> Read %d EVB atoms" % (_MODULE_LABEL, self.natoms))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
