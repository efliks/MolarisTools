#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import TokenizeLine, WriteData
from MolarisAtomsFile    import MolarisAtomsFile

import exceptions, copy, math


_FORMAT_FORCE     = "%16.10f  %16.10f  %16.10f\n"
_FORMAT_CHARGE    = "%16.10f\n"

_FORMAT_SIMPLE    = "%2s   %8.3f   %8.3f   %8.3f\n"
_FORMAT_ARCHIVE   = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n"



class QMCaller (object):
    """Base class to provide communication between Molaris and a QM program.

    This class should not be used directly."""
    # . General options
    # . If qmmm is True, point charges will be used to polarize the wavefunction
    # . Charge schemes can be Mulliken, MerzKollman, Chelpg (only Gaussian)
    defaultAttributes = {
        "method"             :     "PM3"              ,
        "charge"             :      0                 ,
        "multiplicity"       :      1                 ,
        "dielectric"         :     78.4               ,
        "fileAtoms"          :     "mol.in"           ,
        "fileForces"         :     "d.o"              ,
        "fileTrajectory"     :     "qm.xyz"           ,
        "chargeScheme"       :     "Mulliken"         ,
        "archive"            :     False              ,
        "cosmo"              :     False              ,
        "qmmm"               :     False              ,
        "replaceSymbols"     :     [("F0", "F"), ]    ,
        # . Setting this option to True will prevent writing QM forces on MM atoms to the "d.o" file
        "disableQMForces"    :     False              ,
            }

    def __init__ (self, **keywordArguments):
        """Constructor."""
        # . Set default attributes
        attributes = self.__class__.defaultAttributes
        for (key, value) in attributes.iteritems ():
            setattr (self, key, value)

        # . Set user attributes
        for (key, value) in keywordArguments.iteritems ():
            if attributes.has_key (key):
                setattr (self, key, value)
            else:
                raise exceptions.StandardError ("Unknown option %s." % key)

        # . Check for conflicting attributes
        if self.cosmo and self.qmmm:
            raise exceptions.StandardError ("Both cosmo and qmmm options cannot be enabled.")

        # . Preparatory step
        self._Prepare ()


    def _Prepare (self):
        # . Read mol.in file from Molaris
        self.molaris = MolarisAtomsFile (filename=self.fileAtoms, replaceSymbols=self.replaceSymbols)


    def Run (self):
        """Run the calculation."""
        pass


    def _Finalize (self):
        # . Check if the QM calculation is OK
        checks = [ hasattr (self, "Efinal" )  ,
                   hasattr (self, "forces" )  ,
                   hasattr (self, "charges")  ,]
        if not all (checks):
            raise exceptions.StandardError ("Something went wrong.")

        # . Write a file for Molaris containing QM forces and charges
        self._WriteForcesCharges ()
        # . Write a file containing the QM trajectory
        self._WriteTrajectory    ()


    def _WriteForcesCharges (self):
        openfile = open (self.fileForces, "w")
        # . Write the final QM energy
        openfile.write ("%f\n" % self.Efinal)
        # . Write forces on QM atoms
        for force in self.forces:
            openfile.write (_FORMAT_FORCE  % (force.x, force.y, force.z))
        # . Write charges on QM atoms
        for charge in self.charges:
            openfile.write (_FORMAT_CHARGE % charge)
        # . Write forces on MM atoms, if they are present
        if hasattr (self, "mmforces"):
            if not self.disableQMForces:
                for force in self.mmforces:
                    openfile.write (_FORMAT_FORCE  % (force.x, force.y, force.z))
        # . Close the file
        openfile.close ()


    def _WriteTrajectory (self):
        if self.fileTrajectory:
            atoms    = copy.deepcopy (self.molaris.qatoms)
            # . Append link atoms
            atoms.extend (self.molaris.latoms)
            # . Write header
            openfile = open (self.fileTrajectory, "a")
            natoms   = len (atoms)
            openfile.write ("%d\n" % natoms)
            openfile.write ("qm: %f\n" % self.Efinal)

            if self.archive:
                # . Write coordinates, forces and charges
                for atom, force, charge in zip (atoms, self.forces, self.charges):
                    forceMagnitude = math.sqrt (force.x ** 2 + force.y ** 2 + force.z ** 2)
                    openfile.write (_FORMAT_ARCHIVE % (atom.label, atom.x, atom.y, atom.z, force.x, force.y, force.z, forceMagnitude, charge))
            else:
                # . Write coordinates only
                for atom in atoms:
                    openfile.write (_FORMAT_SIMPLE  % (atom.label, atom.x, atom.y, atom.z))
            # . Close the file
            openfile.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
