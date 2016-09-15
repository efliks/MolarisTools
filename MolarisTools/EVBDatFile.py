#-------------------------------------------------------------------------------
# . File      : EVBDatFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine
from    Units     import DEFAULT_EVB_LIB
import  collections, exceptions, os

EVBContainer = collections.namedtuple ("EVBContainer", "serials  types  exist")


_DEFAULT_FILENAME = os.path.join ("evb_heat_01", "evb.dat")
_MODULE_LABEL     = "EVBDat"
_ATOMS_PER_LINE   = 10
_VMD_DIR          = "vmd"


class EVBDatFile (object):
    """A class to represent an evb.dat file."""

    def __init__ (self, filename=_DEFAULT_FILENAME, logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    def _GetContainer (self, data, ntokens):
        serials   = TokenizeLine (next (data), converters=([int, ] * ntokens    ))
        types     = TokenizeLine (next (data), converters=([int, ] * self.nforms))
        exist     = TokenizeLine (next (data), converters=([int, ] * self.nforms))
        container = EVBContainer (
            serials =   serials ,
            types   =   types   ,
            exist   =   exist   , )
        return container


    def _Parse (self, logging):
        data = open (self.filename, "r")
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
        try:
            while True:
                line = next (data)
                if line.count ("# of evb atoms, # of resforms"):
                    (natoms, nforms) = TokenizeLine (line, converters=[int, int])
                    nlines  = (natoms / _ATOMS_PER_LINE) + (1 if (natoms % _ATOMS_PER_LINE) > 0 else 0)
                    # . Read atom serial numbers
                    serials = []
                    for i in range (nlines):
                        tokens = TokenizeLine (next (data), converters=([int, ] * _ATOMS_PER_LINE))
                        for token in tokens:
                            if token != None:
                                serials.append (token)
                    # . Read atom type numbers (?) for each form
                    for i in range (nforms):
                        for j in range (nlines):
                            next (data)
                    # . Read atomic charges for each form
                    charges = []
                    for i in range (nforms):
                        form = []
                        for j in range (nlines):
                            tokens = TokenizeLine (next (data), converters=([float, ] * _ATOMS_PER_LINE))
                            for token in tokens:
                                if token != None:
                                    form.append (token)
                        charges.append (form)
                    self.charges = charges
                    self.natoms  = natoms
                    self.nforms  = nforms
                    if logging:
                        print ("# . %s> Found %d EVB atoms" % (_MODULE_LABEL, natoms))

                elif line.count ("bonds(atoms,types,exist)"):
                    (nbonds, foo) = TokenizeLine (line, converters=[int, None])
                    bonds = []
                    for i in range (nbonds):
                        bond = self._GetContainer (data, 2)
                        bonds.append (bond)
                    self.bonds = bonds
                    if logging:
                        print ("# . %s> Found %d EVB bonds" % (_MODULE_LABEL, nbonds))

                elif line.count ("angles(atoms,types,exist)"):
                    (nangles, foo) = TokenizeLine (line, converters=[int, None])
                    angles = []
                    for i in range (nangles):
                        angle = self._GetContainer (data, 3)
                        angles.append (angle)
                    self.angles = angles
                    if logging:
                        print ("# . %s> Found %d EVB angles" % (_MODULE_LABEL, nangles))

                elif line.count (" torsions(atoms,types,exist)"):
                    (ntorsions, foo) = TokenizeLine (line, converters=[int, None])
                    torsions = []
                    for i in range (ntorsions):
                        torsion = self._GetContainer (data, 4)
                        torsions.append (torsion)
                    self.torsions = torsions
                    if logging:
                        print ("# . %s> Found %d EVB torsions" % (_MODULE_LABEL, ntorsions))
        except StopIteration:
            pass
        # . File closing
        data.close ()


    def GenerateVMDCommands (self, location=_VMD_DIR):
        """Generate a set of VMD commands for measuring bonds, angles and torsions in a QM/MM trajectory."""
        commands = []
        for (i, bond) in enumerate (self.bonds):
            (indexa, indexb) = map (lambda q: (q - 1), bond.serials)
            commands.append ("label  add    Bonds       0/%4d     0/%4d"        % (indexa, indexb))
            (seriala, serialb) = bond.serials
            commands.append ("label  graph  Bonds       %3d     %s/dist_%d_%d.dat" % (i, location, seriala, serialb))

        for (i, angle) in enumerate (self.angles):
            (indexa, indexb, indexc) = map (lambda q: (q - 1), angle.serials)
            commands.append ("label  add    Angles      0/%4d     0/%4d   0/%4d"   % (indexa, indexb, indexc))
            (seriala, serialb, serialc) = angle.serials
            commands.append ("label  graph  Angles      %3d     %s/angl_%d_%d_%d.dat" % (i, location, seriala, serialb, serialc))

        for (i, torsion) in enumerate (self.torsions):
            (indexa, indexb, indexc, indexd) = map (lambda q: (q - 1), torsion.serials)
            commands.append ("label  add    Dihedrals   0/%4d     0/%4d   0/%4d   0/%4d" % (indexa, indexb, indexc, indexd))
            (seriala, serialb, serialc, seriald) = torsion.serials
            commands.append ("label  graph  Dihedrals   %3d     %s/dihe_%d_%d_%d_%d.dat"    % (i, location, seriala, serialb, serialc, seriald))
        return commands


    def _ParseVMDFile (self, filename):
        lines   = open (filename, "r").readlines ()
        collect = []
        for line in lines:
            (foo, value) = TokenizeLine (line, converters=[float, float])
            collect.append (value)
        average = sum (collect) / len (collect)
        return average


    def ReadVMDFiles (self, location=_VMD_DIR):
        """Read .dat files generated by VMD."""
        if hasattr (self, "bonds"):
            collect = []
            for bond in self.bonds:
                (seriala, serialb) = bond.serials
                filename  = os.path.join (location, "dist_%d_%d.dat" % (seriala, serialb))
                average   = self._ParseVMDFile (filename)
                collect.append (average)
            self.averageBonds = collect

        if hasattr (self, "angles"):
            collect = []
            for angle in self.angles:
                (seriala, serialb, serialc) = angle.serials
                filename  = os.path.join (location, "angl_%d_%d_%d.dat" % (seriala, serialb, serialc))
                average   = self._ParseVMDFile (filename)
                collect.append (average)
            self.averageAngles = collect

        if hasattr (self, "torsions"):
            collect = []
            for torsion in self.torsions:
                (seriala, serialb, serialc, seriald) = torsion.serials
                filename  = os.path.join (location, "dihe_%d_%d_%d_%d.dat" % (seriala, serialb, serialc, seriald))
                average   = self._ParseVMDFile (filename)
                collect.append (average)
            self.averageTorsions = collect


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
