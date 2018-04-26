#-------------------------------------------------------------------------------
# . File      : ParametersLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import  collections, os

from MolarisTools.Utilities  import TokenizeLine


Bond        = collections.namedtuple ("Bond"        , "typea  typeb  k  r0")
HBond       = collections.namedtuple ("HBond"       , "typea  typeb  a  b")
Angle       = collections.namedtuple ("Angle"       , "typea  typeb  typec  k  r0")
Torsion     = collections.namedtuple ("Torsion"     , "typeb  typec  k  periodicity  phase")
Improper    = collections.namedtuple ("Improper"    , "atomType  k")
VanDerWaals = collections.namedtuple ("VanDerWaals" , "atomType  attractive  repulsive  mass")

_MODULE_LABEL     = "ParmLib"
_DEFAULT_PARM_LIB = os.path.join (os.environ["HOME"], "DNA_polymerase", "libs", "parm.lib")


class ParametersLibrary (object):
    """A class to represent parameters of the ENZYMIX force field."""

    def __init__ (self, filename=_DEFAULT_PARM_LIB, logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    @property
    def nparameters (self):
        npar = 0
        labels = ("hbonds", "bonds", "angles", "torsions", "impropers", "vdws", )
        for label in labels:
            if hasattr (self, label):
                data = eval ("self.%s" % label)
                npar += len (data)
        return npar

    def __len__ (self):
        return self.nparameters


    def _Parse (self, logging):
        lines = open (self.filename)
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
        try:
            while True:
                line = next (lines)
                if line.count ("BOND PARAMETERS"):
                        # . Read h-bond parameters
                        #   16                                H-BOND PARAMETERS {2A5,2(1X,F10.2)}
                        #   H2   O0  052899.00  016695.00 
                        #   H2   O1  052899.00  016695.00
                        #   (...) 
                        if line.count ("H-BOND"):
                            tokens   = TokenizeLine (line, converters=[int, ])
                            nhbonds  = tokens[0]
                            self.hbonds = {}
                            for i in range (nhbonds):
                                line   = next (lines)
                                tokens = TokenizeLine (line, converters=[None, None, float, float])
                                (typea, typeb, a, b) = tokens
                                bond = HBond (
                                    typea   =   typea   ,
                                    typeb   =   typeb   ,
                                    a       =   a       ,
                                    b       =   b       , )
                                key = "%s-%s" % (typea, typeb)
                                self.hbonds[key] = bond
                            if logging:
                                print ("# . %s> Read %d h-bonds" % (_MODULE_LABEL, nhbonds))


                        # . Read bond parameters
                        #  180                                    BOND PARAMETERS {2A5,2F8.3}
                        #   P5   P5 400.000   6.000
                        #   P2   C8 375.000   3.920
                        #   (...)
                        else:
                            tokens  = TokenizeLine (line, converters=[int, ])
                            nbonds  = tokens[0]
                            self.bonds = {}
                            for i in range (nbonds):
                                line   = next (lines)
                                tokens = TokenizeLine (line, converters=[None, None, float, float])
                                (typea, typeb, k, r0) = tokens
                                bond = Bond (
                                    typea   =   typea   ,
                                    typeb   =   typeb   ,
                                    k       =   k       ,
                                    r0      =   r0      , )
                                key = "%s-%s" % (typea, typeb)
                                self.bonds[key] = bond
                            if logging:
                                print ("# . %s> Read %d bonds" % (_MODULE_LABEL, nbonds))


                # . Read angle parameters
                # 398                                   ANGLE PARAMETERS {3A5,F8.3,F8.1}
                #   P5   P5   P5  140.00   137.0
                #   C2   C2   R*   60.00   109.5
                #   (...)
                elif line.count ("ANGLE PARAMETERS"):
                    tokens  = TokenizeLine (line, converters=[int, ])
                    nangles = tokens[0]
                    self.angles = {}
                    for i in range (nangles):
                        line   = next (lines)
                        tokens = TokenizeLine (line, converters=[None, None, None, float, float])
                        (typea, typeb, typec, k, r0) = tokens
                        angle = Angle (
                            typea   =   typea   ,
                            typeb   =   typeb   ,
                            typec   =   typec   ,
                            k       =   k       ,
                            r0      =   r0      , )
                        key = "%s-%s-%s" % (typea, typeb, typec)
                        self.angles[key] = angle
                    if logging:
                        print ("# . %s> Read %d angles" % (_MODULE_LABEL, nangles))


                # . Read torsion parameters
                #   84                                TORSION PARAMETERS {2A5,2F8.3,F8.1}
                #   P5   P5   0.000   2.000     0.0
                #   C6   C8   1.500   3.000     0.0
                #   (...)
                elif line.count ("TORSION PARAMETERS"):
                    tokens    = TokenizeLine (line, converters=[int, ])
                    ntorsions = tokens[0]
                    self.torsions = {}
                    for i in range (ntorsions):
                        line        = next (lines)
                        tokens      = TokenizeLine (line, converters=[None, None, float, float, float])
                        (typeb, typec, k, periodicity, phase) = tokens
                        # . For some reason, periodicity is written as a float number in the parm.lib file
                        periodicity = int (periodicity)
                        torsion = Torsion (
                            typeb       =   typeb        ,
                            typec       =   typec        ,
                            k           =   k            ,
                            periodicity =   periodicity  ,
                            phase       =   phase        , )
                        key = "%s-%s" % (typeb, typec)
                        self.torsions[key] = torsion
                    if logging:
                        print ("# . %s> Read %d torsions" % (_MODULE_LABEL, ntorsions))


                # . Read improper torsion parameters
                #    6                                I-TORS PARAMETERS {A5,F8.3}
                #   C3   10.000
                #   N3   40.000
                #   (...)
                elif line.count ("I-TORS PARAMETERS"):
                    tokens     = TokenizeLine (line, converters=[int, ])
                    nimpropers = tokens[0]
                    self.impropers = {}
                    for i in range (nimpropers):
                        line     = next (lines)
                        tokens   = TokenizeLine (line, converters=[None, float])
                        (atomType, k) = tokens
                        improper = Improper (
                            atomType = atomType  ,
                            k        = k         , )
                        self.impropers[atomType] = improper
                    if logging:
                        print ("# . %s> Read %d impropers" % (_MODULE_LABEL, nimpropers))


                # . Read van der Waals and mass paramters
                #  130                                VDW AND MASS PARAMETERS {3x,A2,1X,3F9.3}
                #   C6 01956.000  032.000  012.000
                #   C8 01956.000  032.000  012.000
                #   (...)
                elif line.count ("VDW AND MASS PARAMETERS"):
                    tokens   = TokenizeLine (line, converters=[int, ])
                    nparams  = tokens[0]
                    self.nonbonded = {}
                    for i in range (nparams):
                        line    = next (lines)
                        tokens  = TokenizeLine (line, converters=[None, float, float, float])
                        (atomType, repulsive, attractive, mass) = tokens
                        vdw = VanDerWaals (
                            atomType   = atomType    ,
                            repulsive  = repulsive   ,
                            attractive = attractive  ,
                            mass       = mass        , )
                        self.nonbonded[atomType] = vdw
                    if logging:
                        print ("# . %s> Read %d VDW and mass parameters" % (_MODULE_LABEL, nparams))
        except StopIteration:
            pass
        # . Close the file
        lines.close ()


    def GetVDW (self, atomType):
        if hasattr (self, "nonbonded"):
            if (self.nonbonded.has_key (atomType)): return self.nonbonded[atomType]

            generalType = "%s*" % atomType[0]
            if (self.nonbonded.has_key (generalType)): return self.nonbonded[generalType]
        return None

    def GetBond (self, typea, typeb):
        if hasattr (self, "bonds"):
            key = "%s-%s" % (typea, typeb)
            if (self.bonds.has_key (key)): return self.bonds[key]

            yek = "%s-%s" % (typeb, typea)
            if (self.bonds.has_key (yek)): return self.bonds[yek]
        return None

    def GetAngle (self, typea, typeb, typec):
        if hasattr (self, "angles"):
            key = "%s-%s-%s" % (typea, typeb, typec)
            if (self.angles.has_key (key)): return self.angles[key]

            yek = "%s-%s-%s" % (typec, typeb, typea)
            if (self.angles.has_key (yek)): return self.angles[yek]
        return None

    def GetTorsion (self, typeb, typec):
        if hasattr (self, "torsions"):
            key = "%s-%s" % (typeb, typec)
            if (self.torsions.has_key (key)): return self.torsions[key]

            yek = "%s-%s" % (typec, typeb)
            if (self.torsions.has_key (yek)): return self.torsions[yek]
        return None


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__":
    library = ParametersLibrary ()
