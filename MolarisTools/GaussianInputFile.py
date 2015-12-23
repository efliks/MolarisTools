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

    # . Notes:
    # . Memory is in GB
    # . method includes also a basis, for example B3LYP/6-31G*
    # . restart indicates if the wavefunction should be taken from the previous calculation
    # . By default, skip the calculation of CHELPG charges
    defaultAttributes = {
        "ncpu"            :   1          ,
        "memory"          :   1          ,
        "charge"          :   0          ,
        "multiplicity"    :   1          ,
        "qmmm"            :   False      ,
        "cosmo"           :   False      ,
        "method"          :   "PM3"      ,
        "restart"         :   False      ,
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
        qmmm    = "CHARGE=ANGSTROMS"     if self.qmmm    else ""
        cosmo   = "SCRF=(Solvent=Water)" if self.cosmo   else ""
        chelpg  = "POP=CHELPG"           if self.chelpg  else ""
        restart = "GUESS=READ"           if self.restart else ""
        data.append ("# %s %s %s %s %s FORCE NOSYMM\n\n" % (self.method, chelpg, charges, cosmo, restart))
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
