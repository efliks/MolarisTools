#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import WriteData
from MolarisAtomsFile    import MolarisAtomsFile

import subprocess, exceptions, os, glob, shutil


class QMCaller (object):
    """Base class to provide communication between Molaris and QM programs."""

    # . General options
    # . If qmmm is True, point charges will be used to polarize the wavefunction
    defaultAttributes = {
        "method"             :     "PM3"          ,
        "charge"             :     0              ,
        "multiplicity"       :     1              ,
        "qmmm"               :     False          ,
        "cosmo"              :     False          ,
        "dielectric"         :     78.4           ,
        "fileAtoms"          :     "mol.in"       ,
        "fileForces"         :     "forces.out"   ,
        "fileTrajectory"     :     "qm.xyz"       ,
        "filterSymbols"      :     ["MG", "CL", "BR", "Mg", "Cl", "Br", ] ,
        "archive"            :     False          ,
        "directoryArchive"   :     "archive"      ,
            }

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems ():                 setattr (self, key, value)

        if self.cosmo and self.qmmm:
            raise exceptions.StandardError ("Both cosmo and qmmm options cannot be enabled.")
        # . Preparatory step
        self._Prepare ()

        # . Create a directory for archiving
        if self.archive:
            if not os.path.exists (self.directoryArchive):
                os.makedirs (self.directoryArchive)


    def _Prepare (self):
        # . Read atoms.inp from Molaris
        self.molaris = MolarisAtomsFile (filename=self.fileAtoms, filterSymbols=self.filterSymbols)

        # . Write geometries for the purpose of viewing
        if self.fileTrajectory:
            self.molaris.WriteQM (filename=self.fileTrajectory, link=True, append=True)


    def _Archive (self, filename):
        """Archive the log file."""
        if self.archive:
            files = glob.glob (os.path.join (self.directoryArchive, "*.????????"))
            if files:
                lastFile = files[-1]
                counter  = int (lastFile.split (".")[-1])
            else:
                counter  = 0
            shutil.copy2 (filename, os.path.join (self.directoryArchive, "%s.%08d" % (filename, counter + 1)))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
