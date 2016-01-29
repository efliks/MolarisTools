#-------------------------------------------------------------------------------
# . File      : QMCaller.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from Utilities           import WriteData
from MolarisAtomsFile    import MolarisAtomsFile

import exceptions


class QMCaller (object):
    """Base class to provide communication between Molaris and a QM program.

    This class should not be used directly."""
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
        "_logFile"           :     None           ,
        "fileArchive"        :     "archive.log"  ,
            }

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems ():                 setattr (self, key, value)

        if self.cosmo and self.qmmm:
            raise exceptions.StandardError ("Both cosmo and qmmm options cannot be enabled.")
        # . Preparatory step
        self._Prepare ()


    def _Prepare (self):
        # . Read atoms.inp from Molaris
        self.molaris = MolarisAtomsFile (filename=self.fileAtoms, filterSymbols=self.filterSymbols)

        # . Write geometries for the purpose of viewing
        if self.fileTrajectory:
            self.molaris.WriteQM (filename=self.fileTrajectory, link=True, append=True)


    def _Archive (self):
        """Archive the log file."""
        if self.archive:
            logData    = open (self._logFile, "r")
            output     = open (self.fileArchive, "a")
            count      = 1
            try:
                while True:
                    line = logData.next ()
                    if line.count ("Logfile ends"):
                        count += 1
                    output.write (line)
            except StopIteration:
                pass
            width      = 79
            message    = "Logfile ends: %d" % count
            decoration = (width - (len (message) + 2)) / 2
            banner     = " %s %s %s\n\n" % ("@" * decoration, message, "@" * decoration)
            output.write (banner)
            output.close ()
            logData.close ()


    def Run (self):
        """Run the calculation."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
