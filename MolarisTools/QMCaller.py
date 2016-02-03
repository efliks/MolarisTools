#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import WriteData
from MolarisAtomsFile    import MolarisAtomsFile

import exceptions, math, copy


_DEFAULT_METHOD       =  "PM3"
_DEFAULT_CHARGE       =   0
_DEFAULT_MULTIPLICITY =   1
_DEFAULT_DIELECTRIC   =  78.4
_DEFAULT_FILE_ATOMS   =  "mol.in"
_DEFAULT_FILE_FORCES  =  "d.o"
_DEFAULT_FILE_QM_TRAJ =  "qm.xyz"


class QMCaller (object):
    """Base class to provide communication between Molaris and a QM program.

    This class should not be used directly."""
    # . General options
    # . If qmmm is True, point charges will be used to polarize the wavefunction
    defaultAttributes = {
        "method"             :     _DEFAULT_METHOD        ,
        "charge"             :     _DEFAULT_CHARGE        ,
        "multiplicity"       :     _DEFAULT_MULTIPLICITY  ,
        "dielectric"         :     _DEFAULT_DIELECTRIC    ,
        "fileAtoms"          :     _DEFAULT_FILE_ATOMS    ,
        "fileForces"         :     _DEFAULT_FILE_FORCES   ,
        "fileTrajectory"     :     _DEFAULT_FILE_QM_TRAJ  ,
        "archive"            :     False                  ,
        "cosmo"              :     False                  ,
        "qmmm"               :     False                  ,
        "filterSymbols"      :     ["MG", "CL", "BR", "Mg", "Cl", "Br", ] ,
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
        self.molaris = MolarisAtomsFile (filename=self.fileAtoms, filterSymbols=self.filterSymbols)


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
    
            # . Write atoms
            if self.archive:
                # . Write geometry including forces and charges
                for atom, force, charge in zip (atoms, self.forces, self.charges):
                    magnitude = math.sqrt (force.x**2 + force.y**2 + force.z**2)
                    data.append ("%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.4f\n" % (atom.symbol, atom.x, atom.y, atom.z, force.x, force.y, force.z, magnitude, charge))
            else:
                # . Write simple geometry
                for atom in atoms:
                    data.append ("%2s   %8.3f   %8.3f   %8.3f\n" % (atom.symbol, atom.x, atom.y, atom.z))
            # . Update the trajectory file
            WriteData (data, self.fileTrajectory, append=True)


    def _WriteForces (self):
        data = ["%f\n" % self.Efinal]
        # . Write forces
        for force in self.forces:
            data.append ("%14.6f  %14.6f  %14.6f\n" % (force.x, force.y, force.z))

        # . Write charges
        for charge in self.charges:
            data.append ("%14.6f\n" % charge)

        # . Write to a file
        WriteData (data, self.fileForces, append=False)


    def _Finalize (self):
        # . Check if the QM calculation is OK
        checks = [hasattr (self, "Efinal"), hasattr (self, "forces"), hasattr (self, "charges"), ]
        if not all (checks):
            raise exceptions.StandardError ("Something went wrong.")

        # . Write a d.o file
        self._WriteForces ()

        # . Write a qm.xyz file for viewing
        self._WriteQMTrajectory ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
