#-------------------------------------------------------------------------------
# . File      : EVBLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine
import  collections, exceptions


EVBMorsePair = collections.namedtuple ("EVBMorsePair"  ,  "typea  typeb  morseAB  rab  beta  forceHarmonic  radiusHarmonic")
EVBMorseAtom = collections.namedtuple ("EVBMorseAtom"  ,  "evbType  morseD  radius")
EVBAngle     = collections.namedtuple ("EVBAngle"      ,  "evbType  force  angle0  foo  bar")
EVBTorsion   = collections.namedtuple ("EVBTorsion"    ,  "typea  typeb  force  periodicity  phase")
EVBImproper  = collections.namedtuple ("EVBImproper"   ,  "evbType  force  periodicity  angle0")
EVBvdw       = collections.namedtuple ("EVBvdw"        ,  "evbType  repulsive  attractive")



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
                #   (...)
                if   line.startswith ("morse_type"):
                    self.morse = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, morseD, radius = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        atom = EVBMorseAtom (
                            evbType  = evbType  ,
                            morseD   = morseD   ,
                            radius   = radius   , )
                        self.morse.append (atom)
                        line, comment = self._GetLineWithComment (data)


                # EVBMorsePair = collections.namedtuple ("EVBMorsePair"  ,  "typea  typeb  morseAB  rab  beta  fHarmonic  rHarmonic")
                # m_pair         morse_ab    r_ab      beta     f_harmonic   r0_harmonic
                # 'P0'  'O0'      95.0       1.58       2.0        400.0        1.4
                # 'P0'  'O-'     120.0       1.50       2.0        400.0        1.4
                #   (...)
                elif line.startswith ("m_pair"):
                    self.pairs = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        typea, typeb, morseAB, rab, beta, forceHarmonic, radiusHarmonic = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float, float, float, float, float])
                        pair = EVBMorsePair (
                            typea           =  typea          ,
                            typeb           =  typeb          ,
                            morseAB         =  morseAB        ,
                            rab             =  rab            ,
                            beta            =  beta           ,
                            forceHarmonic   =  forceHarmonic  ,
                            radiusHarmonic  =  radiusHarmonic , )
                        self.pairs.append (pair)
                        line, comment = self._GetLineWithComment (data)


                #                  EVB
                # . Read           /     angle parameters
                #          EVB---EVB
                # angle_type     force      angle0
                # 'H0'            0.00        0.0       0.00    1.00
                #   (...)
                elif line.startswith ("angle_type"):
                    self.angles = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, force, angle0, foo, bar = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float, float, float])
                        angle = EVBAngle (
                            evbType  =  evbType  ,
                            force    =  force    ,
                            angle0   =  angle0   ,
                            foo      =  foo      ,
                            bar      =  bar      , )
                        self.angles.append (angle)
                        line, comment = self._GetLineWithComment (data)


                #                  EVB---EVB
                # . Read           /          torsion parameters
                #          EVB---EVB
                # torsion_type   force      periodicity  phase_angle
                # 'P0'  'O0'      3.0        3.0          0.00
                #   (...)
                elif line.startswith ("torsion_type"):
                    self.torsions = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        typea, typeb, force, periodicity, phase = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float, float, float])
                        torsion = EVBTorsion (
                            typea        =  typea       , 
                            typeb        =  typeb       , 
                            force        =  force       , 
                            periodicity  =  periodicity , 
                            phase        =  phase       , )
                        self.torsions.append (torsion)
                        line, comment = self._GetLineWithComment (data)


                #                EVB
                #                /
                # . Read EVB---EVB      improper torsion parameters
                #                \
                #                EVB
                # itorsion_type  force      periodicity  itorsion_angle0
                # 'N-'           10.0       2.0          180.0
                #   (...)
                elif line.startswith ("itorsion_type"):
                    self.impropers = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, force, periodicity, angle0 = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float, float])
                        improper = EVBImproper (
                            evbType      =  evbType     ,
                            force        =  force       ,
                            periodicity  =  periodicity ,
                            angle0       =  angle0      , )
                        self.impropers.append (improper)
                        line, comment = self._GetLineWithComment (data)


                # . Read EVB...EVB nonbonded parameters
                # vdwevb              vdwa      vdwb
                # 'H0'                5.00      0.00
                #   (...)
                elif line.startswith ("vdwevb"):
                    self.vdwevb = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vdw = EVBvdw (
                            evbType    = evbType     ,
                            repulsive  = repulsive   ,
                            attractive = attractive  , )
                        self.vdwevb.append (vdw)
                        line, comment = self._GetLineWithComment (data)


                # . Read EVB...solvent nonbonded parameters
                # solvdw              vdwa      vdwb
                # 'H0'                5.00      0.00
                #   (...)
                elif line.startswith ("solvdw"):
                    self.vdwsol = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vdw = EVBvdw (
                            evbType    = evbType     ,
                            repulsive  = repulsive   ,
                            attractive = attractive  , )
                        self.vdwsol.append (vdw)
                        line, comment = self._GetLineWithComment (data)
        except StopIteration:
            pass
        # . Close the file
        data.close ()


    def List (self, types=["C0", ]):
        """List all energy terms containing a particular atom type."""
        print ("# . Morse atoms")
        for bond in self.morse:
            if bond.evbType in types:
                print ("%2s    %5.1f    %7.2f" % (bond.evbType, bond.morseD, bond.radius))
        print ("# . Morse pairs")
        for pair in self.pairs:
            if pair.typea in types or pair.typeb in types:
                print ("%2s    %2s    %5.1f    %7.2f    %5.1f    %7.1f    %5.1f" % (pair.typea, pair.typeb, pair.morseAB, pair.rab, pair.beta, pair.forceHarmonic, pair.radiusHarmonic))
        print ("# . Angles")
        for angle in self.angles:
            if angle.evbType in types:
                print ("%2s    %5.1f    %7.2f" % (angle.evbType, angle.force, angle.angle0))
        print ("# . Torsions")
        for torsion in self.torsions:
            if torsion.typea in types or torsion.typeb in types:
                print ("%2s    %2s    %7.1f    %7.1f    %7.1f" % (torsion.typea, torsion.typeb, torsion.force, torsion.periodicity, torsion.phase))
        print ("# . Impropers")
        for improper in self.impropers:
            if improper.evbType in types:
                print ("%2s    %7.1f    %7.1f    %7.1f" % (improper.evbType, improper.force, improper.periodicity, improper.angle0))
        print ("# . vdw EVB...EVB")
        for vdw in self.vdwevb:
            if vdw.evbType in types:
                print ("%2s    %10.3f    %10.3f" % (vdw.evbType, vdw.repulsive, vdw.attractive))
        print ("# . vdw EVB...solvent")
        for vdw in self.vdwsol:
            if vdw.evbType in types:
                print ("%2s    %10.3f    %10.3f" % (vdw.evbType, vdw.repulsive, vdw.attractive))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass