#-------------------------------------------------------------------------------
# . File      : MopacInputFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import collections

from MolarisTools.Utilities  import TokenizeLine


Atom = collections.namedtuple ("Atom" , "label x y z")

class MopacInputFile (object):
    """A class to read a MOPAC input file."""

    def __init__ (self, filename="run.mop"):
        """Constructor."""
        self.inputfile = filename
        self._Parse ()


    def _Parse (self):
        lines = open (self.inputfile)
        atoms = []
        try:
            # . Read header
            #   MNDOD  1SCF  CHARGE=-2      GRAD  XYZ  MULLIK  DEBUG MOL_QMMM
            line     = next (lines)
            tokens   = TokenizeLine (line)
            self.keywords = tokens
            xyz      = False
            for token in tokens:
                if   token.startswith ("MNDOD"):
                    pass
                elif token.startswith ("CHARGE"):
                    keyword, charge = token.split ("=")
                    self._charge = float (charge)
                elif token.startswith ("MULLIK"):
                    pass
                elif token.startswith ("XYZ"):
                    xyz = True
            # . Read comment
            line   = next (lines)
            self.comment = line
            # . Skip one line
            line   = next (lines)
            # . Read geometry
            #    C       3.6661  1       6.4850  1      12.9735  1
            #    O       4.8715  1       6.4228  1      13.6710  1
            #    O       4.0062  1       5.5052  1      15.8884  1
            # (...)
            while True:
                line      = next (lines)
                if xyz:
                    tokens    = TokenizeLine (line, converters=[None, float, int, float, int, float, int])
                    # . There may be a blank line at the end
                    if len (tokens) > 0:
                        atomLabel = tokens[0]
                        atomX, atomY, atomZ = tokens[1], tokens[3], tokens[5]
                        atom = Atom (
                            label   =   atomLabel   ,
                            x       =   atomX       ,
                            y       =   atomY       ,
                            z       =   atomZ       ,
                            )
                        atoms.append (atom)
        except StopIteration:
            pass
        # . Close the file
        lines.close ()
        # . Finish up
        self.atoms = atoms


    def Write (self, filename="new.mop"):
        """Write a MOPAC input file."""
        if "XYZ" not in self.keywords:
            raise exceptions.StandardError ("XYZ geometry not defined in the header.")

        fo = open (filename, "w")
        # . Write header, comment, blank line
        header = " ".join (self.keywords) + "\n"
        fo.write (header       )
        fo.write (self.comment )
        fo.write ("\n"         )
        # . Write atoms
        for atom in self.atoms:
            fo.write (" %2s    %9.4f  1    %9.4f  1    %9.4f  1\n" % (atom.label, atom.x, atom.y, atom.z))
        # . Close the file
        fo.close ()


    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0

    @property
    def charge (self):
        if hasattr (self, "_charge"):
            return self._charge
        return 0.


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
