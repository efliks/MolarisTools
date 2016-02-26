#-------------------------------------------------------------------------------
# . File      : AminoLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
# . Atoms are not children of a group.
# . The groups only contain serial numbers of atoms.
#
#               AminoLibrary
#                 /       \
#                /         \
#   AminoComponent     AminoComponent
#                      /      |     \
#                     /       |      \
#              AminoGroups  bonds  AminoAtoms

from    Utilities import TokenizeLine
import  collections, exceptions


AminoAtom  = collections.namedtuple ("Atom"  , "atomNumber  atomLabel  atomType  atomCharge")
AminoGroup = collections.namedtuple ("Group" , "natoms  centralAtom  radius  serials  symbol")

_DEFAULT_DIVIDER = "-" * 41
_GROUP_START     = "A"


class AminoComponent (object):
    """A class to represent a residue."""

    def __init__ (self, logging=True, verbose=False, **keywordArguments):
        """Constructor."""
        for (key, value) in keywordArguments.iteritems ():
            if key != "logging":
                setattr (self, key, value)
        # . Print info
        if logging: self.Info (verbose=verbose)


    def Info (self, verbose=False):
        """Print info."""
        print ("Component: %d %s (%d atoms, %d bonds, %d groups, %-5.2f charge%s)" % (self.serial, self.name, self.natoms, self.nbonds, self.ngroups, self.charge, (", %s" % self.title) if verbose else ""))


    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        else:
            return 0

    @property
    def nbonds (self):
        if hasattr (self, "bonds"):
            return len (self.bonds)
        else:
            return 0

    @property
    def ngroups (self):
        if hasattr (self, "groups"):
            return len (self.groups)
        else:
            return 0

    @property
    def charge (self):
        if hasattr (self, "atoms"):
            total = 0.
            for atom in self.atoms:
                total += atom.atomCharge
            return total
        else:
            return 0.


    def _NumberToLabel (self, atomNumber):
        label = ""
        for atom in self.atoms:
            if atom.atomNumber == atomNumber:
                label = atom.atomLabel
                break
        return label


    def CalculateGroup (self, group):
        """Calculate the charge of a group."""
        total = 0.
        for atom in self.atoms:
            if atom.atomNumber in group.serials:
                total += atom.atomCharge
        return total


    _DEFAULT_MODIFY = {
        "A" :   ("O3",  "O-"),
        "B" :   ("O3",  "O-"),
        "C" :   ("O3",  "O-"),
        }
    def WriteDict (self, groups=None, nstates=2, modify=_DEFAULT_MODIFY):
        """Generate Python code to use in evb_assign.py type of script."""
        # . If no symbols defined, take the whole component
        groupSymbols = groups
        if groupSymbols is None:
            groupSymbols = []
            for group in self.groups:
                groupSymbols.append (group.symbol)
        # . Start here
        print ("names = {")
        for symbol in groupSymbols:
            # . Pick a group of atoms
            found = False
            for group in self.groups:
                if group.symbol == symbol:
                    found = True
                    break
            if not found:
                raise exceptions.StandardError ("Group %s not found." % symbol)
            # . Run for each atom in the group
            for serial in group.serials:
                atom    = self.atoms[serial - 1]
                atype   = atom.atomType
                # . Check if the atom type has to be modified
                if modify.has_key (symbol):
                    oldType, newType = modify[symbol]
                    if oldType == atype:
                        atype = newType
                    else:
                        atype = "%s0" % atype[0]
                else:
                    atype = "%s0" % atype[0]

                acharge = atom.atomCharge
                line    = "%-6s  :  (" % ("\"%s\"" % atom.atomLabel)
                # . For each atom, generate entries for n states
                for i in range (nstates):
                    if i > 0:
                        line = "%s  ,  (%5.2f , \"%2s\")" % (line, acharge, atype)
                    else:
                        line = "%s(%5.2f , \"%2s\")" % (line, acharge, atype)
                print ("%s) ," % line)
            # . Separate groups
            print ("\\")
        # . Finish up
        print ("}")


    def Write (self, title=None, showGroups=False, showLabels=False):
        # . Write header
        if title is None:
            if hasattr (self, "title"):
                title = self.title
        print ("%d%s%s" % (self.serial, self.name, ("  ! %s" % title) if title else ""))

        # . Write atoms
        print ("%5d  ! Number of atoms" % self.natoms)
        for atom in self.atoms:
            markGroup = ""
            if showGroups:
                for igroup, group in enumerate (self.groups):
                    if atom.atomNumber in group.serials:
                        markGroup   = "  %s ! %s" % ("   " if (igroup % 2 == 1) else "", group.symbol)
                        break
            print ("%5d %-4s %4s %6.2f%s" % (atom.atomNumber, atom.atomLabel, atom.atomType, atom.atomCharge, markGroup))

        # . Write bonds
        print ("%5d  ! Number of bonds" % self.nbonds)
        for atoma, atomb in self.bonds:
            label = ""
            if showLabels:
                label = "%4s %4s" % (self._NumberToLabel (atoma), self._NumberToLabel (atomb))
            print ("%5d%5d%s" % (atoma, atomb, ("    ! %s" % label) if showLabels else ""))

        # . Write connecting atoms
        connecta, connectb = self.connect
        print ("%5d%5d  ! Connecting atoms" % (connecta, connectb))

        # . Write groups
        print ("%5d  ! Number of electroneutral groups" % self.ngroups)
        for group in self.groups:
            print ("%5d%5d%6.1f" % (group.natoms, group.centralAtom, group.radius))
            line = "    "
            for atomSerial in group.serials:
                line = "%s%d  " % (line, atomSerial)
            if showGroups:
                line = "%s  ! Group %s: %.4f" % (line, group.symbol, self.CalculateGroup (group))
            print line
        # . Finish up
        print ("%5d" % 0)
        print (_DEFAULT_DIVIDER)


    def KillAtom (self, label):
        """Delete an atom from the component."""
        pass


    def KillBond (self, label, labelOther):
        """Delete a bond from the component."""
        pass


#===============================================================================
class AminoLibrary (object):
    """A class to represent data from the Molaris amino98.lib file."""

    def __init__ (self, filename="amino98_custom.lib", logging=True, reorder=True, unique=False, verbose=False):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging, reorder=reorder, unique=unique, verbose=verbose)


    def _FindComponent (self, key):
        if isinstance (key, int):
            # . Search by serial
            for component in self.components:
                if component.serial == key:
                    return component
        elif isinstance (key, str):
            # . Search by name
            for component in self.components:
                if component.name == key:
                    return component
        else:
            raise exceptions.StandardError ("Unknown type of key.")
        # . Component not found
        return None


    def __contains__ (self, key):
        """Check if a site is in the library."""
        component = self._FindComponent (key)
        if component:
            return True
        return False


    def __getitem__ (self, key):
        """Find and return a component from the library."""
        component = self._FindComponent (key)
        if not component:
            raise exceptions.StandardError ("Site %s not found in the library." % key)
        return component


    @property
    def ncomponents (self):
        if hasattr (self, "components"):
            return len (self.components)
        else:
            return 0


    def __len__ (self):
        return self.ncomponents


    def _GetCleanLine (self, data):
        line      = data.next ()
        lineClean = line[:line.find ("!")].strip ()
        return lineClean


    def _GetLineWithComment (self, data):
        line      = data.next ()
        position  = line.find ("!")
        text      = line[             : position].strip ()
        comment   = line[position + 1 :         ].strip ()
        return (text, comment)


    def _Parse (self, logging, reorder, unique, verbose=False):
        components = []
        names      = []
        data       = open (self.filename)
        try:
            while True:
                line = self._GetCleanLine (data)
                # . Check if a new residue starts
                if line.startswith ("---"):
                    # . Get serial and name
                    line, title  = self._GetLineWithComment (data)
                    entry = TokenizeLine (line, converters=[None, ])[0]
                    # . Remove spaces
                    entry = entry.replace (" ", "")
                    # . Check if last residue found
                    if entry == "0":
                        break
                    for i, char in enumerate (entry):
                        if not char.isdigit ():
                            break
                    serial, name = int (entry[:i]), entry[i:]
                    # . Check if the component label is unique
                    if unique:
                        if name in names:
                            raise exceptions.StandardError ("Component label %s is not unique." % name)
                        names.append (name)
                    # . Get number of atoms
                    line   = self._GetCleanLine (data)
                    natoms = int (line)
                    # . Read atoms
                    atoms  = []
                    labels = []
                    for i in range (natoms):
                        line = self._GetCleanLine (data)
                        atomNumber, atomLabel, atomType, atomCharge = TokenizeLine (line, converters=[int, None, None, float])
                        if unique:
                            # . Check if the atom label is unique
                            if atomLabel in labels:
                                raise exceptions.StandardError ("Component %s %d: Atom label %s is not unique." % (name, serial, atomLabel))
                            labels.append (atomLabel)
                        # . Create atom
                        atom = AminoAtom (atomNumber=atomNumber, atomLabel=atomLabel, atomType=atomType, atomCharge=atomCharge)
                        atoms.append (atom)
                    # . Get number of bonds
                    line   = self._GetCleanLine (data)
                    nbonds = int (line)
                    # . Read bonds
                    bonds  = []
                    for i in range (nbonds):
                        line = self._GetCleanLine (data)
                        atoma, atomb = TokenizeLine (line, converters=[int, int])
                        if reorder:
                            # . Keep the lower number first
                            if atoma > atomb:
                                atoma, atomb = atomb, atoma
                        bonds.append ((atoma, atomb))
                    if reorder:
                        # . Sort bonds
                        bonds.sort (key=lambda bond: bond[0])
                    # . Read connecting atoms
                    line = self._GetCleanLine (data)
                    connecta, connectb = TokenizeLine (line, converters=[int, int])
                    # . Read number of electroneutral groups
                    line    = self._GetCleanLine (data)
                    ngroups = int (line)
                    # . Read groups
                    groups  = []
                    for i in range (ngroups):
                        line     = self._GetCleanLine (data)
                        nat, central, radius = TokenizeLine (line, converters=[int, int, float])
                        line     = self._GetCleanLine (data)
                        serials  = TokenizeLine (line, converters=[int] * nat)
                        symbol   = chr (ord (_GROUP_START) + i)
                        group    = AminoGroup (natoms=nat, centralAtom=central, radius=radius, serials=serials, symbol=symbol)
                        groups.append (group)
                    # . Create a component and add it to the list
                    component = AminoComponent (serial=serial, name=name, atoms=atoms, bonds=bonds, groups=groups, connect=(connecta, connectb), logging=logging, title=title, verbose=verbose)
                    components.append (component)
        except StopIteration:
            pass
        # . Finish up
        data.close ()
        if logging:
            print ("Found %d components." % len (components))
        self.components = components


    def WriteComponent (self, serial=None, name=None, title="", showGroups=False, showLabels=False):
        """Write out a specific component from the library."""
        component = self.GetComponent (serial=serial, name=name)
        component.Write (title=title, showGroups=showGroups, showLabels=showLabels)


    def WriteAll (self, showGroups=False, showLabels=False):
        """Write out all components from the library."""
        for component in self.components:
            component.Write (showGroups=showGroups, showLabels=showLabels)
        # . Write footer
        print ("%5d" % 0)
        print (_DEFAULT_DIVIDER)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
