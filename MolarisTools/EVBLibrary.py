#-------------------------------------------------------------------------------
# . File      : EVBLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine
import  collections, exceptions


class EVBLibrary (object):
    """A class to represent a collection of EVB paramters."""

    def __init__ (self, filename="evb_poll_clean.lib", logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    def _GetLineWithComment (self, data):
        line      = data.next ()
        position  = line.find ("!")
        if position > -1:
            text      = line[             : position]
            comment   = line[position + 1 :         ]
        else:
            text      = line
            comment   = ""
        text    = text.strip ()
        comment = comment.strip ()
        return (text, comment)


    def _Parse (self, logging):
        data = open (self.filename)
        try:
            while True:
                line, comment = self._GetLineWithComment (data)
                # . Read EVB---EVB bond parameters
                # morse_type     morse_d    radius
                # 'H0'           105.0      0.4
                # ...
                if   line.startswith ("morse_type"):
                    self.bonds = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        evbType, morse_d, radius = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        bond = (evbType, morse_d, radius)
                        self.bonds.append (bond)


                #                  EVB
                # . Read           /     angle parameters
                #          EVB---EVB
                # angle_type     force      angle0
                # 'H0'            0.00        0.0       0.00    1.00
                # ...
                elif line.startswith ("angle_type"):
                    self.angles = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        evbType, force, angle0, foo, bar = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float, float, float])
                        angle = (evbType, force, angle0, foo, bar)
                        self.angles.append (angle)


                #                  EVB---EVB
                # . Read           /          torsion parameters
                #          EVB---EVB
                # torsion_type   force      periodicity  phase_angle
                # 'P0'  'O0'      3.0        3.0          0.00
                # ...
                elif line.startswith ("torsion_type"):
                    self.torsions = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        typea, typeb, force, periodicity, phase = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float, float, float])
                        torsion = (typea, typeb, force, periodicity, phase)
                        self.torsions.append (torsion)


                #                EVB
                #                /
                # . Read EVB---EVB      improper torsion parameters
                #                \
                #                EVB
                # itorsion_type  force      periodicity  itorsion_angle0
                # 'N-'           10.0       2.0          180.0
                # ...
                elif line.startswith ("itorsion_type"):
                    self.impropers = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        evbType, force, periodicity, angle = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float, float])
                        improper = (evbType, force, periodicity, angle)
                        self.impropers.append (improper)


                # . Read EVB...EVB nonbonded parameters
                # vdwevb              vdwa      vdwb
                # 'H0'                5.00      0.00
                # ...
                elif line.startswith ("vdwevb"):
                    self.vdwevb = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vdw = (evbType, repulsive, attractive)
                        self.vdwevb.append (vdw)


                # . Read EVB...solvent nonbonded parameters
                # solvdw              vdwa      vdwb
                # 'H0'                5.00      0.00
                # ...
                elif line.startswith ("solvdw"):
                    self.vdwsol = []
                    while True:
                        line, comment = self._GetLineWithComment (data)
                        if line == "":
                            break
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vdw = (evbType, repulsive, attractive)
                        self.vdwsol.append (vdw)
        except StopIteration:
            pass
        # . Close the file
        data.close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
