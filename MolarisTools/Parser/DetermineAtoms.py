#-------------------------------------------------------------------------------
# . File      : DetermineAtoms.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import  collections, exceptions

from MolarisTools.Library  import AminoLibrary


Residue = collections.namedtuple ("Residue" , "label  serial  charge  atoms")
Atom    = collections.namedtuple ("Atom"    , "label  serial  charge  connectLabels  connectSerials  x  y  z")

class DetermineAtoms (object):
    """A class to read a Molaris output file containing residue information."""

    def __init__ (self, filename="determine_atoms.out"):
        """Constructor."""
        self.filename = filename
        self._Parse ()


    def _Parse (self):
        """Parse output file."""
        data     = open (self.filename)
        residues = []
        try:
            while True:
                line = next (data)
                if line.count ("atom list for residue"):
                    # . Found a new residue
                    #
                    # atom list for residue:   247_ASP,    # of atoms in this residue:  12
                    tokens    = line.split ()
                    residue   = tokens[4]
                    resSerial, resLabel = residue.split ("_")
                    resSerial = int (resSerial)
                    resLabel  = resLabel.replace (",", "")
                    # . Skip a few lines
                    for i in range (3):
                        next (data)


                    # . Read atom information
                    #
                    # number  name  type      x          y          z      charge      atoms bonded(name)      atoms bonded(number)
                    # ------  ----  ----   -------    -------    -------   ------   ------------------------ ------------------------
                    # 3935    N      N3      0.048     10.092     19.698   -0.400   C    HN   CA              3933  3936  3937
                    # 3936    HN     H2     -0.378      9.269     19.322    0.400   N                         3935
                    # 3937    CA     C4      0.489     11.137     18.785   -0.097   N    HA   CB   C          3935  3938  3939  3945
                    # ...
                    atoms  = []
                    natoms = int (tokens[-1])
                    for i in range (natoms):
                        line         = next (data)
                        tokens       = line.split ()
                        (atomSerial, atomLabel, atomType), atomCharge = tokens[:3], tokens[6]
                        atomSerial   = int (atomSerial)
                        atomCharge   = float (atomCharge)
                        items        = tokens[7:]
                        bondLabels   = []
                        bondSerials  = []
                        for item in items:
                            if item.isdigit ():
                                serial = int (item)
                                bondSerials.append (serial)
                            else:
                                bondLabels.append (item)
                        # . Get coordinates (useful for imposing positional restraints)
                        x, y, z = map (float, tokens[3:6])
                        # . Create a new atom
                        atom = Atom (label=atomLabel, serial=atomSerial, charge=atomCharge, connectLabels=bondLabels, connectSerials=bondSerials, x=x, y=y, z=z)
                        atoms.append (atom)
                    # . Skip one line
                    next (data)


                    # . Extract the total charge calculated by Molaris
                    #
                    # Total charge of this residue:    -1.000
                    line        = next (data)
                    tokens      = line.split ()
                    totalCharge = float (tokens[-1])

                    # . Create a residue
                    residue = Residue (label=resLabel, serial=resSerial, charge=totalCharge, atoms=atoms)
                    residues.append (residue)
        except StopIteration:
            pass
        # . Close the file
        data.close ()
        # . Finish up
        self.residues = residues


    @property
    def nresidues (self):
        if hasattr (self, "residues"):
            return len (self.residues)
        else:
            return 0


    _DEFAULT_MODIFY = {
            "O3" : "O-" ,
            "MG" : "MG" ,
        }
    def GenerateList (self, selectResidues, library, nstates=2, shift=2, modify=_DEFAULT_MODIFY):
        """Generate a list of EVB atoms to use in a Molaris input file."""
        if not isinstance (library, AminoLibrary):
            raise exceptions.StandardError ("Argument is not an amino library.")

        # . Iterate residues
        residues = []
        for selectLabel, selectSerial, selectGroups in selectResidues:
            # . Perform initial checks
            if not selectLabel in library:
                raise exceptions.StandardError ("Residue %s not found in the library." % selectLabel)

            # . Iterate residues from the determine_atoms.out file
            found = False
            for selfResidue in self.residues:
                checks = (
                    selfResidue.label  == selectLabel  ,
                    selfResidue.serial == selectSerial ,
                         )
                if all (checks):
                    found = True
                    break
            if not found:
                raise exceptions.StandardError ("Residue %s not found in file %s." % (selectLabel, self.filename))


            # . Select AminoComponent
            aminoComponent = library[selectLabel]

            # . Select atoms belonging to selected groups
            selectedAtoms  = []
            for aminoGroup in aminoComponent.groups:
                if aminoGroup.symbol in selectGroups:
                    selectedAtoms.extend (aminoGroup.labels)
            aminoAtoms     = []
            for aminoAtom in aminoComponent.atoms:
                if aminoAtom.atomLabel in selectedAtoms:
                    aminoAtoms.append (aminoAtom)


            # . Print a list of EVB atoms
#        evb_atm      7196    0.7000    P0    0.7000    P0      #     0.7000  P4    PG  
#        evb_atm      7197   -0.9000    O-   -0.9000    O-      #    -0.9000  O3    O1G 
#        evb_atm      7198   -0.9000    O-   -0.9000    O-      #    -0.9000  O3    O2G 
            for selfAtom in selfResidue.atoms:
                for aminoAtom in aminoAtoms:
                    if selfAtom.label == aminoAtom.atomLabel:
                        # . Modify atom type
                        if modify.has_key (aminoAtom.atomType):
                            modifiedType = modify[aminoAtom.atomType]
                        else:
                            modifiedType = "%s0" % aminoAtom.atomType[0]
                        entry  = "     %7.4f    %2s" % (aminoAtom.atomCharge, modifiedType)
                        states = ""
                        for i in range (nstates):
                            states = "%s%s" % (states, entry)
                        print ("%sevb_atm     %5d%s       #  %7.4f   %2s     %-4s" % ("    " * shift, selfAtom.serial, states, aminoAtom.atomCharge, aminoAtom.atomType, aminoAtom.atomLabel))
                        break

            # . Print a list of EVB bonds
#        evb_bnd   0      7203  7200   # C3B     PB    
#        evb_bnd   0      7201  7200   # O1B     PB    
#        evb_bnd   0      7202  7200   # O2B     PB    


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
