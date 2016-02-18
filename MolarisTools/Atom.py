#-------------------------------------------------------------------------------
# . File      : Atom.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import  math

_FORMAT_QM_SIMPLE        = "%2s   %8.3f   %8.3f   %8.3f\n"
_FORMAT_QM_FORCE         = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n"
_FORMAT_QM_FORCE_CHARGE  = "%2s   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.4f\n"


class Atom (object):
    """A class to handle an atom."""

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, value) in keywordArguments.iteritems ():
            setattr (self, key, value)
        # label  x  y  z  fx  fy  fz  fm  charge

        # . Calculate the magnitude of the force
        checks = (
            not hasattr (self, "fm") ,
                hasattr (self, "fx") ,)
        if all (checks):
            self.fm = math.sqrt (self.fx ** 2 + self.fy ** 2 + self.fz ** 2)

    @property
    def line (self):
        if hasattr (self, "fx"):
            if hasattr (self, "charge"):
                return _FORMAT_QM_FORCE_CHARGE % (self.label, self.x, self.y, self.z, self.fx, self.fy, self.fz, self.fm, self.charge)
            else:
                return _FORMAT_QM_FORCE % (self.label, self.x, self.y, self.z, self.fx, self.fy, self.fz, self.fm)
        else:
            return _FORMAT_QM_SIMPLE % (self.label, self.x, self.y, self.z)
