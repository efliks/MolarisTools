#-------------------------------------------------------------------------------
# . File      : AminoLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Utilities import TokenizeLine
import collections, exceptions

Atom  = collections.namedtuple ("Atom"  , "atomNumber  atomLabel  atomType  atomCharge")
Group = collections.namedtuple ("Group" , "natoms  centralAtom  radius  atoms")

DEFAULT_DIVIDER = "-" * 41

# TODO
# 1. Checking if charges in a group sum up to integer values, write to stderr if they don't
#


class AminoComponent (object):
    """A class to represent a residue."""

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in keywordArguments.iteritems (): setattr (self, key, value)

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


    def Write (self, title="", showGroups=False, showLabels=False):
        # . Write header
        print ("%d%s%s" % (self.serial, self.name, ("  ! %s" % title) if title else ""))

        # . Write atoms
        print ("%5d  ! Number of atoms" % self.natoms)
        for atom in self.atoms:
            markGroup = ""
            if showGroups:
                for igroup, group in enumerate (self.groups):
                    if atom.atomNumber in group.atoms:
                        groupSymbol = chr (ord ("A") + igroup)
                        markGroup   = "    %s ! %s" % (" " * igroup * 2, groupSymbol)
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
            for atomSerial in group.atoms:
                line = "%s%d  " % (line, atomSerial)
            print line
        # . Finish up
        print ("%5d" % 0)
        print (DEFAULT_DIVIDER)


#===============================================================================
class AminoLibrary (object):
    """A class to represent data from the Molaris amino98.lib file."""

    def __init__ (self, filename="amino98_custom.lib", logging=True, reorder=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging, reorder=reorder)


    @property
    def ncomponents (self):
        if hasattr (self, "components"):
            return len (self.components)
        else:
            return 0


    def _GetCleanLine (self, data):
        line      = data.next ()
        lineClean = line[:line.find ("!")].strip ()
        return lineClean


    def _Parse (self, logging, reorder):
        components = []
        data       = open (self.filename)
        try:
            while True:
                line = self._GetCleanLine (data)
                # . Check if a new residue starts
                if line.startswith ("---"):
                    # . Get serial and name
                    line  = self._GetCleanLine (data)
                    entry = TokenizeLine (line, converters=[None, ])[0]
                    # . Check if last residue found
                    if entry == "0":
                        break
                    for i, char in enumerate (entry):
                        if not char.isdigit ():
                            break
                    serial, name = int (entry[:i]), entry[i:]
                    # . Get number of atoms
                    line   = self._GetCleanLine (data)
                    natoms = int (line)
                    # . Read atoms
                    atoms  = []
                    for i in range (natoms):
                        line = self._GetCleanLine (data)
                        atomNumber, atomLabel, atomType, atomCharge = TokenizeLine (line, converters=[int, None, None, float])
                        atom = Atom (atomNumber=atomNumber, atomLabel=atomLabel, atomType=atomType, atomCharge=atomCharge)
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
                        elements = TokenizeLine (line, converters=[int] * nat)
                        group    = Group (natoms=nat, centralAtom=central, radius=radius, atoms=elements)
                        groups.append (group)
                    if logging:
                        print ("Found component: %d %s (%d atoms, %d bonds, %d groups)" % (serial, name, natoms, nbonds, ngroups))
                    # . Create a component and add it to the list
                    component = AminoComponent (serial=serial, name=name, atoms=atoms, bonds=bonds, groups=groups, connect=(connecta, connectb))
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
        found = False
        for component in self.components:
            if component.serial == serial or component.name == name:
                found = True
                break
        if found:
            component.Write (title=title, showGroups=showGroups, showLabels=showLabels)
        else:
            raise exceptions.StandardError ("Component not found.")


    def WriteAll (self, showGroups=False, showLabels=False):
        """Write out all components from the library."""
        for component in self.components:
            component.Write (showGroups=showGroups, showLabels=showLabels)
        # . Write footer
        print ("%5d" % 0)
        print (DEFAULT_DIVIDER)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
