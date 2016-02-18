#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import TokenizeLine, WriteData
from Atom                import Atom
from XYZTrajectory       import XYZStep
from MolarisAtomsFile    import MolarisAtomsFile

import exceptions, copy


_FORMAT_FORCE     = "%14.6f  %14.6f  %14.6f\n"
_FORMAT_CHARGE    = "%14.6f\n"


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
            atoms    = copy.deepcopy (self.molaris.qatoms)
            atoms.extend (self.molaris.latoms)

            if self.archive:
                templ = []
                for atom, force, charge in zip (atoms, self.forces, self.charges):
                    atom = Atom (
                        label  = atom.label ,
                        charge = charge     ,
                        x      = atom.x     ,
                        y      = atom.y     ,
                        z      = atom.z     ,
                        fx     = force.x    ,
                        fy     = force.y    ,
                        fz     = force.z    ,)
                    templ.append (atom)
                atoms = templ
            # . Write a step
            step = XYZStep (atoms, comment="qm: %f" % self.Efinal)
            step.Write (filename=self.fileTrajectory, append=True)


    def _WriteForces (self):
        data = ["%f\n" % self.Efinal]
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
