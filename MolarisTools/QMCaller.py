#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import WriteData
from MolarisAtomsFile    import MolarisAtomsFile

import exceptions, math, copy


_DEFAULT_METHOD         =  "PM3"
_DEFAULT_CHARGE         =   0
_DEFAULT_MULTIPLICITY   =   1
_DEFAULT_DIELECTRIC     =  78.4
_DEFAULT_FILE_ATOMS     =  "mol.in"
_DEFAULT_FILE_FORCES    =  "d.o"
_DEFAULT_FILE_QM_TRAJ   =  "qm.xyz"
# . Other schemes: MerzKollman, Chelpg (only Gaussian)
_DEFAULT_CHARGE_SCHEME  =  "Mulliken"

_FORMAT_QM_EXT    = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.4f\n"
_FORMAT_QM_SIMPLE = "%2s   %8.3f   %8.3f   %8.3f\n"
_FORMAT_FORCE     = "%14.6f  %14.6f  %14.6f\n"
_FORMAT_CHARGE    = "%14.6f\n"


class QMCaller (object):
    """Base class to provide communication between Molaris and a QM program.

    This class should not be used directly."""
    # . General options
    # . If qmmm is True, point charges will be used to polarize the wavefunction
    defaultAttributes = {
        "method"           :   _DEFAULT_METHOD         ,
        "charge"           :   _DEFAULT_CHARGE         ,
        "multiplicity"     :   _DEFAULT_MULTIPLICITY   ,
        "dielectric"       :   _DEFAULT_DIELECTRIC     ,
        "fileAtoms"        :   _DEFAULT_FILE_ATOMS     ,
        "fileForces"       :   _DEFAULT_FILE_FORCES    ,
        "fileTrajectory"   :   _DEFAULT_FILE_QM_TRAJ   ,
        "chargeScheme"     :   _DEFAULT_CHARGE_SCHEME  ,
        "archive"          :   False                   ,
        "cosmo"            :   False                   ,
        "qmmm"             :   False                   ,
        "replaceSymbols"   :   [("F0", "F"), ]         ,
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


    def _WriteQMTrajectory (self):
        if self.fileTrajectory:
            # . Collect atoms
            atoms = copy.deepcopy (self.molaris.qatoms)
            atoms.extend (self.molaris.latoms)
    
            # . Write header
            caption = "qm: %f" % self.Efinal
            natoms  = len (atoms)
            data    = ["%d\n%s\n" % (natoms, caption)]
            if self.archive:
                # . Write geometry including forces and charges
                for atom, force, charge in zip (atoms, self.forces, self.charges):
                    magnitude = math.sqrt (force.x * force.x + force.y * force.y + force.z * force.z)
                    data.append (_FORMAT_QM_EXT % (atom.symbol, atom.x, atom.y, atom.z, force.x, force.y, force.z, magnitude, charge))
            else:
                # . Write simple geometry
                for atom in atoms:
                    data.append (_FORMAT_QM_SIMPLE % (atom.symbol, atom.x, atom.y, atom.z))
            # . Update the trajectory file
            WriteData (data, self.fileTrajectory, append=True)


    def _WriteForces (self):
        data = ["%f\n" % self.Efinal]
        # . Write forces and charges
        for force in self.forces:
            data.append (_FORMAT_FORCE % (force.x, force.y, force.z))
        for charge in self.charges:
            data.append (_FORMAT_CHARGE % charge)

        # . Write to a file
        WriteData (data, self.fileForces, append=False)


    def _Finalize (self):
        # . Check if the QM calculation is OK
        checks = [hasattr (self, "Efinal"), hasattr (self, "forces"), hasattr (self, "charges"), ]
        if not all (checks):
            raise exceptions.StandardError ("Something went wrong.")

        # . Write forces for Molaris
        self._WriteForces ()

        # . Write QM trajectory for viewing
        self._WriteQMTrajectory ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
