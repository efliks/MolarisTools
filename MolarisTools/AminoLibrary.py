#-------------------------------------------------------------------------------
# . File      : AminoLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
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
from    Units           import DEFAULT_AMINO_LIB
import  exceptions, os

_MODULE_LABEL    = "AminoLib"
_GROUP_START     = "A"
_DEFAULT_TOPOLOGY_FORMAT = "Molaris"


class AminoLibrary (object):
    """A class to represent data from the Molaris amino98.lib file.

    Can also read CHARMM topology files."""

    def __init__ (self, filename=DEFAULT_AMINO_LIB, logging=True, reorder=True, unique=False, verbose=False, topologyFormat=_DEFAULT_TOPOLOGY_FORMAT, cutType=True):
        """Constructor."""
        self.filename = filename
        if   topologyFormat == "Molaris":
            self._Parse (logging=logging, reorder=reorder, unique=unique, verbose=verbose)
        elif topologyFormat == "CHARMM":
            self._ParseCHARMM (logging=logging, cutType=cutType)
        else:
            raise exceptions.StandardError ("Unknown topology format.")
        if logging:
            print ("# . %s> Found %d component%s" % (_MODULE_LABEL, self.ncomponents, "s" if self.ncomponents != 1 else ""))
            if not verbose:
                self.WriteLabels ()
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


    def _ParseCHARMM (self, logging, cutType):
        data       = open (self.filename)
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
        # . Initialize
        bonds      = []
        group      = []
        groups     = []
        components = []
        try:
            while True:
                line = self._GetCleanLine (data)
                #  RESI ASP         -1.00
                #  GROUP
                #  ATOM N    NH1    -0.47
                #  ATOM H    H       0.31
                #  ATOM CA   CT1     0.07
                #  ATOM HA   HB      0.09
                #  GROUP
                #  ATOM CB   CT2    -0.28
                #   (...)
                #  BOND CB CA  CG CB  OD2 CG
                #   (...)
                #  DOUBLE  O   C   CG  OD1
                #
                # . Treat patches as components
                if line[:4] in ("RESI", "PRES", ):
                    if group:
                        groups.append (group)
                        # . Create a temporary component
                        component = (componentLabel, componentCharge, groups, bonds)
                        components.append (component)
                    # . Component begins
                    tokens    = TokenizeLine (line, converters=[None, None, float])
                    (componentLabel, componentCharge) = tokens[1:]
                    # . Reinitialize
                    bonds     = []
                    group     = []
                    groups    = []
                elif line.startswith ("GROUP"):
                    if group:
                        groups.append (group)
                        group = []
                elif line.startswith ("ATOM" ):
                    tokens  = TokenizeLine (line, converters=[None, None, None, float])
                    newAtom = AminoAtom (
                        # . Molaris only uses 2-character atom types
                        atomType    =    tokens[2][:2] if cutType else tokens[2] ,
                        atomLabel   =    tokens[1]  ,
                        atomCharge  =    tokens[3]  ,
                        )
                    group.append (newAtom)
                elif line.startswith ("BOND" ) or line.startswith ("DOUBLE"):
                    tokens   = line.split ()
                    labels   = tokens[1:]
                    nlabels  = len (labels)
                    if (nlabels % 2) != 0:
                        raise exceptions.StandardError ("Incorrect BOND entry in component %s." % componentLabel)
                    for i in range (0, nlabels, 2):
                        labela, labelb = labels[i], labels[i + 1]
                        # . Ignore bonds involving atoms from other residues
                        checks = []
                        for label in (labela, labelb):
                            checks.extend ( [label[ 0] != "+", label[-1] != "-", not label[ 0].isdigit ()] )
                        if all (checks):
                            pair = (labela, labelb)
                            bonds.append (pair)
        except StopIteration:
            pass
        # . Finish up
        data.close ()
        if group:
            groups.append (group)
            component = (componentLabel, componentCharge, groups, bonds)
            components.append (component)
        # . Set up actual amino components from temporary components
        aminoComponents = []
        for componentSerial, (componentLabel, componentCharge, groups, bonds) in enumerate (components, 1):
            # . Merge atoms
            aminoAtoms  = []
            for group in groups:
                for atom in group:
                    aminoAtoms.append (atom)
            aminoGroups = []
            # . Iterate temporary groups
            for igroup, group in enumerate (groups):
                # . Merge atom labels
                labels = []
                for atom in group:
                    labels.append (atom.atomLabel)
                # . Create a group
                natoms = len (group) 
                aminoGroup = AminoGroup (
                    radius      =   5.       ,
                    natoms      =   natoms   ,
                    labels      =   labels   ,
                    symbol      =   chr (ord (_GROUP_START) + igroup) ,
                    centralAtom =   group[natoms / 2].atomLabel       ,
                    )
                aminoGroups.append (aminoGroup)
            # . Create a component
            component = AminoComponent (
                serial  =   componentSerial ,
                label   =   componentLabel  ,
                groups  =   aminoGroups     ,
                atoms   =   aminoAtoms      ,
                bonds   =   bonds           ,
                connect =   ("",    "")     ,
                logging =   logging         ,
                title   =   "Generated from CHARMM topology"  ,
                )
            aminoComponents.append (component)
        self.components = aminoComponents


    def _Parse (self, logging, reorder, unique, verbose):
        components = []
        names      = []
        data       = open (self.filename)
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
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
                    logFlag = False
                    if logging:
                        if verbose:
                            logFlag = True
                    component = AminoComponent (serial=componentSerial, name=name, atoms=atoms, bonds=bonds, groups=groups, connect=(connecta, connectb), logging=logFlag, title=title)
                    components.append (component)
        except StopIteration:
            pass
        # . Finish up
        data.close ()
        self.components = components


    def WriteAll (self, showGroups=False, showLabels=False):
        """Write out all components from the library."""
        for component in self.components[:-1]:
            component.Write (showGroups=showGroups, showLabels=showLabels, terminate=False)

        # . Write the last component
        component = self.components[-1]
        component.Write (showGroups=showGroups, showLabels=showLabels, terminate=True)


    def WriteLabels (self, rowl=14, serials=False):
        """Write labels of all components."""
        line = ""
        for (i, component) in enumerate (self.components, 1):
            if serials:
                line = "%s   %3d %-3s" % (line, component.serial, component.label)
            else:
                line = "%s  %3s" % (line, component.label)
            if (i % rowl == 0):
                print line
                line = ""
        if line:
            print line


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__":
    filename = os.path.join (os.environ["HOME"], "devel", "pcetk", "tests", "lysozyme", "toppar", "top_all27_prot_na.inp")
    library  = AminoLibrary (filename=filename, logging=True, topologyFormat="CHARMM")

    component = library["ARG"]
    component.Write (showGroups=True, showLabels=True, sortGroups=True, terminate=True)
