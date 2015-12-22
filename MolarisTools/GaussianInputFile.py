#-------------------------------------------------------------------------------
# . File      : GaussianInputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from   Utilities    import WriteData
import collections

Atom  = collections.namedtuple ("Atom"  , "symbol x y z charge")


class GaussianInputFile (object):
    """A class for writing Gaussian input files."""

    defaultAttributes = {
        "ncpu"            :   1          ,
        # . memory in GB
        "memory"          :   1          ,
        "charge"          :   0          ,
        "multiplicity"    :   1          ,
        "qmmm"            :   False      ,
        # . method includes also a basis, for example B3LYP/6-31G*
        "method"          :   "PM3"      ,
        # . restart indicates if the wavefunction should be taken from the previous calculation
        "restart"         :   False      ,
        # . By default, skip the calculation of CHELPG charges
        "chelpg"          :   False      ,
        "fileInput"       :   "job.inp"  ,
        "fileCheckpoint"  :   None       ,
            }

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems ():                 setattr (self, key, value)


    def Write (self, atoms, pointCharges=None):
        """Write a Gaussian input file from a Molaris atoms file."""
        data = []
        # . Write job control
        if self.ncpu > 1:
            data.append ("%%NProcShared=%d\n" % self.ncpu)
        if self.memory:
            data.append ("%%mem=%dgb\n" % self.memory)
        if self.fileCheckpoint:
            data.append ("%%chk=%s\n" % self.fileCheckpoint)

        # . Write header
        charges = "CHARGE=ANGSTROMS" if self.qmmm    else ""
        restart = "GUESS=READ"       if self.restart else ""
        chelpg  = "POP=CHELPG"       if self.chelpg  else ""
        data.append ("# %s %s FORCE %s %s NOSYMM\n\n" % (self.method, chelpg, charges, restart))
        data.append ("Comment line\n\n")
        data.append ("%d %d\n" % (self.charge, self.multiplicity))

        # . Write geometry
        for atom in atoms:
            data.append ("%2s    %9.4f    %9.4f    %9.4f\n" % (atom.symbol, atom.x, atom.y, atom.z))
        data.append ("\n")

        # . Write point charges
        if pointCharges:
            if self.qmmm:
                for atom in pointCharges:
                    data.append ("%9.4f    %9.4f    %9.4f    %9.4f\n" % (atom.x, atom.y, atom.z, atom.charge))
                data.append ("\n")

        # . Finish up
        WriteData (data, self.fileInput)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
