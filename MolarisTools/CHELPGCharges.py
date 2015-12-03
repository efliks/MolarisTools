#-------------------------------------------------------------------------------
# . File      : CHELPGCharges.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
"""Extract and manipulate CHELPG charges from a Gaussian or ORCA log file."""

import os
import sys
import exceptions

class Atom (object):
#    def __init__ (self, symbol, atomicNumber, x, y, z, charge, chargeOrig):
#        self.symbol       = symbol
#        self.atomicNumber = atomicNumber
#        self.charge       = charge
#        self.chargeOrig   = chargeOrig
#        self.x            = x
#        self.y            = y
#        self.z            = z
# http://stackoverflow.com/questions/2280334/shortest-way-of-creating-an-object-with-arbitrary-attributes-in-python
    def __init__(self, **kw):
        for name in kw:
            setattr (self, name, kw[name])


class CHELPGCharges (object):
    """A class for reading CHELPG charges from output files."""

    def __init__ (self, logFile):
        """Constructor."""
        self.title   = os.path.basename (logFile)
        self.logFile = logFile
        if self._DetectFormat () == "ORCA":
            self._ParseORCA ()
        else:
            self._ParseGaussian ()


    def _DetectFormat (self):
        lines = open (self.logFile).readlines ()
        frm   = "Gaussian"
        for line in lines:
            if line.count ("O   R   C   A"):
                frm = "ORCA"
                break
        return frm


    def _ParseGaussian (self):
        """A parser of Gaussian files."""
        lines      = open (self.logFile).readlines ()
        reading    = False
        # . Read in the charges
        for line in lines:
            if not reading:
                if line.startswith (" Charges from ESP fit"):
                    atoms   = []
                    reading = True
            else:
                if line.startswith (" ---"):
                    reading = False
                else:
                    tokens = line.split ()
                    tokn   = len (tokens)
                    if tokn == 3:
                        serial, symbol, charge = tokens
                        charge = float (charge)
                        atoms.append ([serial, symbol, charge])
        # . Read in the geometry
        reading = False
        for line in lines:
            if not reading:
                if line.count ("Standard orientation"):
                    geom    = []
                    count   = 0
                    reading = True
            else:
                count += 1
                if line.startswith (" ---") and count > 4:
                    reading = False
                else:
                    tokens = line.split ()
                    if tokens[0].isdigit ():
                        serial, atomicNumber, atomicType, x, y, z = tokens
                        x, y, z = map (float, [x, y, z])
                        geom.append ([atomicNumber, x, y, z])
        # . Merge arrays
        self.atoms  = []
        for (serial, symbol, charge), (number, x, y, z) in zip (atoms, geom):
            atom = Atom (symbol=symbol, atomicNumber=number, x=x, y=y, z=z, charge=charge, chargeOrig=charge)
            self.atoms.append (atom)


    def _ParseORCA (self):
        """A parser of ORCA files."""
        lines      = open (self.logFile).readlines ()
        reading    = False
        itr        = iter (lines)
        try:
            while True:
                line = itr.next ()
                # . Read in the charges
                if line.startswith ("CHELPG Charges"):
                    atoms = []
                    itr.next ()
                    while True:
                        line   = itr.next ()
                        if line.startswith ("---"):
                            break
                        tokens = line.split ()
                        index, symbol, foo, charge = tokens
                        index  = int (index)
                        charge = float (charge)
                        atoms.append ([index + 1, symbol, charge])
                # . Read in the geometry
                elif line.startswith ("CARTESIAN COORDINATES (ANGSTROEM)"):
                    geom = []
                    itr.next ()
                    while True:
                        line   = itr.next ()
                        tokens = line.split ()
                        if len (tokens) < 2:
                            break
                        symbol, x, y, z = tokens
                        x, y, z = map (float, [x, y, z])
                        geom.append ([-1, x, y, z])
        except StopIteration:
            pass
        # . Merge arrays
        self.atoms  = []
        for (serial, symbol, charge), (number, x, y, z) in zip (atoms, geom):
            atom = Atom (symbol=symbol, atomicNumber=number, x=x, y=y, z=z, charge=charge, chargeOrig=charge)
            self.atoms.append (atom)


    def AverageCharges (self, serials):
        """Average some charges, for example for hydrogen atoms in the methyl group."""
        # . Collect individual charges
        total = 0.
        for serial in serials:
            atom   = self.atoms[serial - 1]
            total += atom.charge
        # . Replace individual charges by the average charge
        average = total / len (serials)
        for serial in serials:
            atom        = self.atoms[serial - 1]
            atom.charge = average


    def MergeCharges (self, serials):
        """Set the charges of some atoms to zero and add the removed charge to the charge of another atom."""
        # . Collect individual charges and zet them to zero
        total = 0.
        for serial in serials[:-1]:
            atom         = self.atoms[serial - 1]
            total       += atom.charge
            atom.charge  = 0.
        # . Merge the collected charges with the charge of the last atom
        lastSerial   = serials[-1]
        atom         = self.atoms[lastSerial - 1]
        atom.charge += total


    def GroupCharges (self, serials, groupCharge, serialsOther=None):
        """Create a group of charges.

        Currently, it is possible to define only one group per molecule."""
        # . Calculate the current total charge of the group
        total    = 0.
        nserials = len (serials)
        for serial in serials:
            atom   = self.atoms[serial - 1]
            total += atom.charge
        # . Update the charges of the group atoms
        delta        = groupCharge - total
        deltaPerAtom = delta / nserials
        for serial in serials:
            atom         = self.atoms[serial - 1]
            atom.charge += deltaPerAtom
        # . Update the charges of other atoms so that the total charge of the molecule remains the same
        natoms       = len (self.atoms) - nserials
        deltaPerAtom = -delta / natoms
        for serial in range (1, len (self.atoms) + 1):
            if serial not in serials:
                atom         = self.atoms[serial - 1]
                atom.charge += deltaPerAtom
        # . Remember the group
        self.group = serials


    def FixCharge (self, serial, charge):
        """Set a predefined charge to an atom.

        Charges have to be fixed before other operations are performed.

        Currently, fixed charges and groups of charges cannot be used at the same time."""
        # . Maintain a list of fixed atoms
        if hasattr (self, "fixed"):
            self.fixed.append (serial)
        else:
            self.fixed = [serial]
        # . Set the new charge
        atom         = self.atoms[serial - 1]
        oldCharge    = atom.charge
        atom.charge  = charge
        # . Update charges of the other atoms
        delta        = charge - oldCharge
        deltaPerAtom = -delta / (len (self.atoms) - len (self.fixed))
        for otherSerial, otherAtom in enumerate (self.atoms, 1):
            if otherSerial not in self.fixed:
                otherAtom.charge += deltaPerAtom


    def WriteGeometryCharges (self, extraCol=False, extraLine=False, line=""):
        """Write an XYZ file with charges."""
        total     = 0.
        totalOrig = 0.
        for atom in self.atoms:
            total     += atom.charge
            totalOrig += atom.chargeOrig
        natoms = len (self.atoms)

        print ("%4d" % natoms)
        print ("%s %.2f %.2f" % (self.title, total, totalOrig))
        if extraCol:
            for serial, atom in enumerate (self.atoms, 1):
                group = ""
                if hasattr (self, "group"):
                    group = "*" if serial in self.group else " "
                if hasattr (self, "fixed"):
                    group = "*" if serial in self.fixed else " "
                print ("%2s  %8.3f %8.3f %8.3f    %6.2f    # %1s%-2d   %6.2f" % (atom.symbol, atom.x, atom.y, atom.z, atom.charge, group, serial, atom.chargeOrig))
        else:
            for atom in self.atoms:
                print ("%2s  %8.3f %8.3f %8.3f    %6.2f" % (atom.symbol, atom.x, atom.y, atom.z, atom.charge))
        if extraLine:
            for atom in self.atoms:
                line = "%s %6.2f" % (line, atom.charge)
            print line
