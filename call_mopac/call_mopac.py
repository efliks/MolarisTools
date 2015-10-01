#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : call_mopac.py
# . Copyright : USC, Mikolaj J. Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections, exceptions, subprocess

Atom  = collections.namedtuple ("Atom"  , "symbol x y z charge")
Force = collections.namedtuple ("Force" , "x y z")

# http://users.mccammon.ucsd.edu/~dzhang/energy-unit-conv-table.html
HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM = 627.509469 / 0.529177
HARTREE_TO_KCAL_MOL               = 627.509469
GRADIENT_TO_FORCE                 =  -1.
EV_TO_KCAL_MOL                    =  23.0609


def _TokenizeLine (line, converters=None, separator=None, reverse=False):
    """Tokenize a line with optional converters and separator."""
    tokens = None
    if line is not None:
        if separator is None:
            tokens = line.split ()
        else:
            tokens = line.split (separator)
        # . The last token becomes first
        if reverse:
            tokens.reverse ()
        if converters is not None:
            # . Truncate tokens if there is fewer of them than converters
            ntokens     = len (tokens)
            nconverters = len (converters)
            if ntokens > nconverters:
                tokens = tokens[:nconverters]
            for (i, (token, converter)) in enumerate (zip (tokens, converters)):
                if converter is None:
                    new = token
                else:
                    try:
                        new = converter (token)
                    except:
                        raise exceptions.StandardError ("Unable to convert token " + repr (i) + ".", True)
                tokens[i] = new
    return tokens

def WriteData (data, filename, append=False):
        output = open (filename, "a" if append else "w")
        output.writelines (data)
        output.close ()


#-------------------------------------------------------------------------------
class MolarisFile (object):
    """A class representing atoms for the QC/MM calculation."""

    def __init__ (self, filename="atoms.inp"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _LineToAtom (self, line, includeCharge=False):
        if includeCharge:
            tokens = _TokenizeLine (line, converters=[None, float, float, float, float])
            atom   = Atom (symbol=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], charge=tokens[4])
        else:
            tokens = _TokenizeLine (line, converters=[None, float, float, float])
            atom   = Atom (symbol=tokens[0], x=tokens[1], y=tokens[2], z=tokens[3], charge=None)
        return atom


    def _ReadAtoms (self, openfile, natoms, includeCharge=False):
        atoms = []
        for nq in range (natoms):
            line = next (openfile)
            atom = self._LineToAtom (line, includeCharge)
            atoms.append (atom)
        return atoms


    def _Parse (self):
        lines   = open (self.inputfile)
        line    = next (lines)
        # . Get step number and total energy
        step, G = _TokenizeLine (line, converters=[int, float])
        try:
            while True:
                line = next (lines)

                # . Read the QM section
                if line.count ("# of qmmm atoms"):
                    nquantum, nlink = _TokenizeLine (line, converters=[int, int])
                    # . Read QM atoms proper
                    self.qatoms = self._ReadAtoms (lines, nquantum)
                    # . Read QM link atoms
                    self.latoms = self._ReadAtoms (lines, nlink)

                elif line.count ("# of total frozen protein atoms, # of groups in Region I`"):
                    pass
                elif line.count ("# of frozen water atoms in Region I`"):
                    pass

                # . Read the protein section
                elif line.count ("# of non-frozen protein atoms in Region II"):
                    nprot = _TokenizeLine (line, converters=[int])[0]
                    self.patoms = self._ReadAtoms (lines, nprot, includeCharge=True)

                # . Read the free water section
                elif line.count ("# of non-frozen water atoms in the system"):
                    nwater = _TokenizeLine (line, converters=[int])[0]
                    self.watoms = self._ReadAtoms (lines, nwater, includeCharge=True)
        except StopIteration:
            pass


    def WriteQM (self, filename="qm.xyz", link=False, append=False):
        """Write QM atoms to an XYZ file."""
        self._WriteSection (filename, (self.qatoms + self.latoms) if link else self.qatoms, includeCharge=False, append=append)


    def WriteProtein (self, filename="prot.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.patoms, includeCharge=True)


    def WriteWater (self, filename="wat.xyz"):
        """Write protein atoms to an XYZ file."""
        self._WriteSection (filename, self.watoms, includeCharge=True)


    def _WriteSection (self, filename, atoms, includeCharge=False, skipSymbol=False, append=False):
        """Write an XYZ file."""
        natoms = len (atoms)
        data   = ["%d\n%s\n" % (natoms, filename)]
        for atom in atoms:
            data.append ("%2s    %8.3f    %8.3f    %8.3f    %8s\n" % ("" if skipSymbol else atom.symbol, atom.x, atom.y, atom.z, ("%8.4f" % atom.charge) if includeCharge else ""))
        WriteData (data, filename, append=append)


    def WriteMopacInput (self, filename="run.mop", method="PM3", charge=0):
        """Write an input file for MOPAC."""
        data   = []
        data.append ("%s  1SCF  CHARGE=%-2d  GRAD  XYZ  ESP\n" % (method, charge))
        data.append ("Comment line\n")
        data.append ("\n")
        for atom in self.qatoms + self.latoms:
            data.append ("%2s    %9.4f  1    %9.4f  1    %9.4f  1\n" % (atom.symbol, atom.x, atom.y, atom.z))
        WriteData (data, filename)


    def WriteGaussianInput  (self, filename="run_gauss.inp", methodBasis="PM3", charge=0, multiplicity=1):
        """Write an input file for Gaussian."""
        # . Determine if a semiempirical potential is used
        isSemiempirical = methodBasis[:3] in ("AM1", "PM3", )

        data = []
        data.append ("# %s POP=CHELPG FORCE %s\n\n" % ((methodBasis, "" if isSemiempirical else "CHARGE=ANGSTROMS NOSYMM")))
        data.append ("Comment line\n\n")
        data.append ("%d %d\n" % (charge, multiplicity))
        for atom in self.qatoms + self.latoms:
            data.append ("%2s    %9.4f    %9.4f    %9.4f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        data.append ("\n")

        # . Write point charges
        if not isSemiempirical:
            for atom in (self.patoms + self.watoms):
                data.append ("%9.4f    %9.4f    %9.4f    %9.4f\n" % (atom.x, atom.y, atom.z, atom.charge))
            data.append ("\n")
        WriteData (data, filename)


#-------------------------------------------------------------------------------
class GaussianFile (object):
    """A class to read a Gaussian output file."""

    def __init__ (self, filename="run_gauss.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get the number of atoms
                if line.count ("Input orientation:") or line.count ("Z-Matrix orientation:"):
                    for skip in range (4):
                        next (lines)
                    natoms = 0
                    while True:
                        line = next (lines)
                        if line.count ("----"):
                            break
                        natoms += 1

                # . Get the final energy (for semiempirical calculations)
                elif line.count ("NIter="):
                    tokens      = _TokenizeLine (line, converters=[None, float, None, None])
                    self.Efinal = tokens[1] * HARTREE_TO_KCAL_MOL

                # . Get the final energy (for ab initio/DFT calculations)
                elif line.count ("SCF Done"):
                    tokens      = _TokenizeLine (line, converters=[None, None, None, None, float, None, None, int, None])
                    self.Efinal = tokens[4] * HARTREE_TO_KCAL_MOL

                # . Get ESP charges
                elif line.count ("Charges from ESP fit"):
                    next (lines)
                    next (lines)
                    self.charges = []
                    for i in range (natoms):
                        tokens = _TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.charges.append (charge)

                # . Get Mulliken charges
                elif line.count ("Mulliken atomic charges:"):
                    next (lines)
                    self.mcharges = []
                    for i in range (natoms):
                        tokens = _TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.mcharges.append (charge)

                # . Get forces
                # . http://www.gaussian.com/g_tech/g_ur/k_force.htm
                # elif line.count ("***** Axes restored to original set *****"):
                #     for skip in range (4):
                #         next (lines)
                elif line.count ("Forces (Hartrees/Bohr)"):
                    next (lines)
                    next (lines)
                    self.forces = []
                    for i in range (natoms):
                        tokens = _TokenizeLine (next (lines), converters=[int, int, float, float, float])
                        force  = Force (x=tokens[2] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, y=tokens[3] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM, z=tokens[4] * HARTREE_BOHR_TO_KCAL_MOL_ANGSTROM)
                        self.forces.append (force)
        except StopIteration:
            pass


    def WriteMolarisForces (self, filename="forces.out", Eref=0.):
        """Write a file in the Molaris-suitable format."""
        data = []
        data.append ("%f\n" % (self.Efinal - Eref))
        for force in self.forces:
            data.append ("%14.6f  %14.6f  %14.6f\n" % (force.x, force.y, force.z))
        for charge in self.charges:
            data.append ("%14.6f\n" % charge)

        # . Write the file
        WriteData (data, filename)


#-------------------------------------------------------------------------------
class MopacFile (object):
    """A class to read a MOPAC output file."""

    def __init__ (self, filename="run.out"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _GetGradientLine (self, openfile):
        line     = next (openfile)
        tokens   = _TokenizeLine (line, converters=[int, int, None, None, None, float, float, None])
        gradient = tokens[6]
        return gradient


    def _Parse (self):
        lines = open (self.inputfile)
        try:
            while True:
                line = next (lines)
                # . Get the number of atoms
                if line.count ("TOTAL NO. OF ATOMS:"):
                    tokens = _TokenizeLine (line, converters=[int], reverse=True)
                    natoms = tokens[0]

                # . Read gradients
                elif line.count ("FINAL  POINT  AND  DERIVATIVES"):
                    # . Skip the next two lines
                    next (lines)
                    next (lines)
                    self.forces = []
                    for i in range (natoms):
                        force = Force (x=self._GetGradientLine (lines) * GRADIENT_TO_FORCE, y=self._GetGradientLine (lines) * GRADIENT_TO_FORCE, z=self._GetGradientLine (lines) * GRADIENT_TO_FORCE)
                        self.forces.append (force)

                # . Get the final total energy
                # elif line.count ("TOTAL ENERGY"):
                #     tokens = _TokenizeLine (line, converters=[None, None, None, float, None])
                #     self.Etotal = tokens[3] * EV_TO_KCAL_MOL

                # . Read the final QM energy
                elif line.count ("FINAL HEAT OF FORMATION"):
                    tokens = _TokenizeLine (line, converters=[None, None, None, None, None, float, None, None, float, None])
                    self.Efinal = tokens [5]

                # . Read ESP charges
                elif line.count ("ELECTROSTATIC POTENTIAL CHARGES"):
                    # . Skip the next two lines
                    next (lines)
                    next (lines)
                    self.charges = []
                    for i in range (natoms):
                        tokens = _TokenizeLine (next (lines), converters=[int, None, float])
                        charge = tokens[2]
                        self.charges.append (charge)
        except StopIteration:
            pass


    def WriteMolarisForces (self, filename="forces.out", Eref=0.):
        """Write a file in the Molaris-suitable format."""
        data = []
        data.append ("%f\n" % (self.Efinal - Eref))
        for force in self.forces:
            data.append ("%14.6f  %14.6f  %14.6f\n" % (force.x, force.y, force.z))
        for charge in self.charges:
            data.append ("%14.6f\n" % charge)

        # . Write the file
        WriteData (data, filename)


#===============================================================================
# . Main program
#===============================================================================
# . Read atoms.inp from Molaris
molaris = MolarisFile (filename="atoms.inp")

# . Write geometries for the purpose of viewing
molaris.WriteQM (link=True, append=True)

# . Run a MOPAC calculation
molaris.WriteMopacInput (filename="run.mop", method="PM3", charge=-5)
mopacError  = open ("run.err", "w")
subprocess.check_call (["/home/mikolaj/local/bin/MOPAC2009.exe", "run.mop"], stdout=mopacError, stderr=mopacError)
mopacError.close ()

# . Convert the output file from MOPAC to forces.out
mopac = MopacFile (filename="run.out")
mopac.WriteMolarisForces (filename="forces.out")

# . Run a Gaussian calculation
# molaris.WriteGaussianInput (filename="run_gauss.inp", methodBasis="B3LYP/6-31G(d)", charge=-5)

# ... Run Gaussian here ...

# . Convert the output file from Gaussian to forces.out
# gaussian = GaussianFile (filename="run_gauss.out")
# gaussian.WriteMolarisForces (filename="forces_gauss.out")
