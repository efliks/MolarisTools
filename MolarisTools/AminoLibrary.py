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

from    Utilities       import TokenizeLine
from    AminoComponent  import AminoAtom, AminoGroup, AminoComponent
import  exceptions, os

_MODULE_LABEL    = "AminoLib"
_GROUP_START     = "A"
_DEFAULT_LIBRARY_FILE   = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "amino98_custom_small.lib")


class AminoLibrary (object):
    """A class to represent data from the Molaris amino98.lib file."""

    def __init__ (self, filename=_DEFAULT_LIBRARY_FILE, logging=True, reorder=True, unique=False, verbose=False):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging, reorder=reorder, unique=unique, verbose=verbose)
        self._i = 0


    def __len__ (self):
        return self.ncomponents


    # . The next 3 methods are for the iterator
    def __iter__ (self):
        return self

    def __next__ (self):
        return self.next ()

    def next (self):
        """Next component."""
        if self._i >= self.ncomponents:
            self._i = 0
            raise exceptions.StopIteration ()
        else:
            self._i += 1
        return self.components[self._i - 1]


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


    def has_key (self, key):
        """Checks for a component in the library."""
        return self.__contains__ (key)


    def __contains__ (self, key):
        """Checks for a component in the library."""
        component = self._FindComponent (key)
        if component:
            return True
        return False


    def __getitem__ (self, key):
        """Find and return a component from the library."""
        component = self._FindComponent (key)
        if not component:
            raise exceptions.StandardError ("Component %s not found in the library." % key)
        return component


    @property
    def ncomponents (self):
        if hasattr (self, "components"):
            return len (self.components)
        else:
            return 0


    @property
    def lastSerial (self):
        serial = 1
        if self.ncomponents > 1:
            for component in self.components:
                serial = component.serial
        return serial


    def _GetCleanLine (self, data):
        line      = data.next ()
        lineClean = line[:line.find ("!")].strip ()
        return lineClean


    def _GetLineWithComment (self, data):
        line      = data.next ()
        position  = line.find ("!")
        if position > -1:
            text      = line[             : position].strip ()
            comment   = line[position + 1 :         ].strip ()
        else:
            text      = line
            comment   = ""
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
                    componentSerial, name = int (entry[:i]), entry[i:]
                    # . Check if the component label is unique
                    if unique:
                        if name in names:
                            raise exceptions.StandardError ("Component label %s is not unique." % name)
                        names.append (name)
                    # . Get number of atoms
                    line    = self._GetCleanLine (data)
                    natoms  = int (line)
                    # . Initiate conversion table serial->label
                    convert = {}
                    # . Read atoms
                    atoms   = []
                    labels  = []
                    for i in range (natoms):
                        line = self._GetCleanLine (data)
                        atomNumber, atomLabel, atomType, atomCharge = TokenizeLine (line, converters=[int, None, None, float])
                        if unique:
                            # . Check if the atom label is unique
                            if atomLabel in labels:
                                raise exceptions.StandardError ("Component %s %d: Atom label %s is not unique." % (name, serial, atomLabel))
                            labels.append (atomLabel)
                        # . Create atom
                        atom = AminoAtom (atomLabel=atomLabel, atomType=atomType, atomCharge=atomCharge)
                        atoms.append (atom)
                        # . Update conversion table serial->label
                        convert[atomNumber] = atomLabel
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
                    # . Convert numerical bonds to labeled bonds
                    labeledBonds = []
                    for atoma, atomb in bonds:
                        # . FIXME: Workaround for invalid entries in the amino file
                        try:
                            pair = (convert[atoma], convert[atomb])
                            labeledBonds.append (pair)
                        except:
                            pass
                    bonds = labeledBonds
                    # . Read connecting atoms
                    line  = self._GetCleanLine (data)
                    seriala, serialb = TokenizeLine (line, converters=[int, int])
                    # . Convert serials of connecting atoms to labels
                    connecta, connectb = "", ""
                    if seriala > 0:
                        connecta = convert[seriala]
                    if serialb > 0:
                        connectb = convert[serialb]
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
                        # . Convert central atom's serial to a label
                        # . FIXME: Workaround for invalid entries in the amino file
                        try:
                            central  = convert[central]
                            # . Convert serials to labels
                            labels   = []
                            for serial in serials:
                                labels.append (convert[serial])
                            symbol   = chr (ord (_GROUP_START) + i)
                            group    = AminoGroup (natoms=nat, centralAtom=central, radius=radius, labels=labels, symbol=symbol)
                            groups.append (group)
                        except:
                            pass
                    # . Create a component and add it to the list
                    component = AminoComponent (serial=componentSerial, name=name, atoms=atoms, bonds=bonds, groups=groups, connect=(connecta, connectb), logging=logging, title=title, verbose=verbose)
                    components.append (component)
        except StopIteration:
            pass
        # . Finish up
        data.close ()
        if logging:
            ncomponents = len (components)
            print ("# . %s> Found %d component%s" % (_MODULE_LABEL, ncomponents, "s" if ncomponents > 1 else ""))
        self.components = components


    def WriteAll (self, showGroups=False, showLabels=False):
        """Write out all components from the library."""
        for component in self.components[:-1]:
            component.Write (showGroups=showGroups, showLabels=showLabels, terminate=False)

        # . Write the last component
        component = self.components[-1]
        component.Write (showGroups=showGroups, showLabels=showLabels, terminate=True)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
