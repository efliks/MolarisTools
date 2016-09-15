#-------------------------------------------------------------------------------
# . File      : EVBLibrary.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2016)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities import TokenizeLine
from    Units     import DEFAULT_EVB_LIB
import  collections, exceptions, math, os


EVBMorsePair  = collections.namedtuple ("EVBMorsePair"  ,  "typea  typeb  morseAB  rab  beta  forceHarmonic  radiusHarmonic")
EVBMorseAtom  = collections.namedtuple ("EVBMorseAtom"  ,  "evbType  morseD  radius")
EVBAngle      = collections.namedtuple ("EVBAngle"      ,  "evbType  force  angle0  foo  bar")
EVBTorsion    = collections.namedtuple ("EVBTorsion"    ,  "typea  typeb  force  periodicity  phase")
EVBImproper   = collections.namedtuple ("EVBImproper"   ,  "evbType  force  periodicity  angle0")
EVBvdw        = collections.namedtuple ("EVBvdw"        ,  "evbType  repulsive  attractive")
EVBNonBond    = collections.namedtuple ("EVBNonBond"    ,  "evbType  exRepulsive  beta  qvdwa  qvdwb")
EVBPair       = collections.namedtuple ("EVBPair"       ,  "typea  typeb  exRepulsive  beta  qvdwa  qvdwb")
EVBInductive  = collections.namedtuple ("EVBInductive"  ,  "evbType  alpha  screen")
EVBIndPair    = collections.namedtuple ("EVBIndPair"    ,  "typea  typeb  alpha")
EVBScreen     = collections.namedtuple ("EVBScreen"     ,  "evbType  mus")
EVBScreenPair = collections.namedtuple ("EVBScreenPair" ,  "typea  typeb  mus")

_MODULE_LABEL     = "EVBLib"
_COMMENT_CHARS    = ("!", "#")


class EVBLibrary (object):
    """A class to represent a collection of EVB paramters."""

    def __init__ (self, filename=DEFAULT_EVB_LIB, logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    def _GetLineWithComment (self, data):
        line      = next (data)
        positions = []
        for char in _COMMENT_CHARS:
            position = line.find (char)
            if position >= 0:
                positions.append (position)
        if positions:
            positions.sort ()
            position  = positions[0]
            text      = line[             : position]
            comment   = line[position + 1 :         ]
        else:
            text      = line
            comment   = ""
        text    = text.strip ()
        comment = comment.strip ()
        return (text, comment)


    def _Parse (self, logging):
        messages = {
            "morse_type"        :   "# . %s> Read %d Morse types"   ,
            "m_pair"            :   "# . %s> Read %d Morse pairs"   ,
            "angle_type"        :   "# . %s> Read %d angles"        ,
            "torsion_type"      :   "# . %s> Read %d torsions"      ,
            "itorsion_type"     :   "# . %s> Read %d impropers"     ,
            "vdwevb"            :   "# . %s> Read %d EVB...EVB VDW parameters"      ,
            "solvdw"            :   "# . %s> Read %d EVB...solvent VDW parameters"  ,
            "exnonbond_type"    :   "# . %s> Read %d nonbonded EVB...EVB parameters for atom types"     ,
            "x_pair"            :   "# . %s> Read %d nonbonded EVB...EVB parameters for atom pairs"     ,
            "v_pair"            :   "# . %s> Read %d nonbonded EVB...EVB parameters for atom types (never bonded)"  ,
            "induct"            :   "# . %s> Read %d inductive interactions by atom type"               ,
            "a_induct"          :   "# . %s> Read %d inductive interactions by atom pair"               ,
            "elect"             :   "# . %s> Read %d electrostatic screening corrections by atom type"  ,
            "a_elect"           :   "# . %s> Read %d electrostatic screening corrections by atom pair"  , }
        data = open (self.filename)
        npar = 0
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
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
                    if logging:
                        nmorse = len (self.morse)
                        npar += nmorse
                        print (messages["morse_type"] % (_MODULE_LABEL, nmorse))


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
                    if logging:
                        npairs = len (self.pairs)
                        npar += npairs
                        print (messages["m_pair"] % (_MODULE_LABEL, npairs))


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
                    if logging:
                        nangles = len (self.angles)
                        npar += nangles
                        print (messages["angle_type"] % (_MODULE_LABEL, nangles))


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
                    if logging:
                        ntorsions = len (self.torsions)
                        npar += ntorsions
                        print (messages["torsion_type"] % (_MODULE_LABEL, ntorsions))


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
                    if logging:
                        nimpropers = len (self.impropers)
                        npar += nimpropers
                        print (messages["itorsion_type"] % (_MODULE_LABEL, nimpropers))


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
                    if logging:
                        nvdw = len (self.vdwevb)
                        npar += nvdw
                        print (messages["vdwevb"] % (_MODULE_LABEL, nvdw))


                # . Read EVB...solvent nonbonded parameters
                # solvdw              vdwa      vdwb
                # 'H0'                5.00      0.00
                #   (...)
                elif line.startswith ("solvdw"):
                    self.solvdw = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vdw = EVBvdw (
                            evbType    = evbType     ,
                            repulsive  = repulsive   ,
                            attractive = attractive  , )
                        self.solvdw.append (vdw)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nvdw = len (self.solvdw)
                        npar += nvdw
                        print (messages["solvdw"] % (_MODULE_LABEL, nvdw))


                # . Nonbonded EVB...EVB interaction (by atom type)
                # exnonbond_type  ex_repl   beta     q_vdwa    q_vdwb
                # 'H0'            5.00      2.50 	  0.00000     0.000
                # (...)
                elif line.startswith ("exnonbond_type"):
                    self.exnonbond = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, exRepulsive, beta, qvdwa, qvdwb = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float, float, float])
                        exnonbond = EVBNonBond (
                            evbType     = evbType      ,
                            exRepulsive = exRepulsive  ,
                            beta        = beta         ,
                            qvdwa       = qvdwa        ,
                            qvdwb       = qvdwb        , )
                        self.exnonbond.append (exnonbond)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nexnb = len (self.exnonbond)
                        npar += nexnb
                        print (messages["exnonbond_type"] % (_MODULE_LABEL, nexnb))


                # . Nonbonded EVB...EVB interaction (by atom pair)
                # x_pair          ex_repl    beta     q_vdwa    q_vdwb
                # 'L0' 'C0'        69999.0      4.00     0.     0.
                # (...)
                elif line.startswith ("x_pair"):
                    self.xpairs = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        typea, typeb, exRepulsive, beta, qvdwa, qvdwb = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float, float, float, float])
                        xpair = EVBPair (
                            typea       = typea        ,
                            typeb       = typeb        ,
                            exRepulsive = exRepulsive  ,
                            beta        = beta         ,
                            qvdwa       = qvdwa        ,
                            qvdwb       = qvdwb        , )
                        self.xpairs.append (xpair)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nxpair = len (self.xpairs)
                        npar += nxpair
                        print (messages["x_pair"] % (_MODULE_LABEL, nxpair))


                # . Nonbonded EVB...EVB interaction for never bonded atoms
                # v_pair              vdwa      vdwb
                # (...)
                elif line.startswith ("v_pair"):
                    self.vpairs = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, repulsive, attractive = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        vpair = EVBvdw (
                            evbType    = evbType     ,
                            repulsive  = repulsive   ,
                            attractive = attractive  , )
                        self.vpairs.append (vpair)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nvpair = len (self.vpairs)
                        npar += nvpair
                        print (messages["v_pair"] % (_MODULE_LABEL, nvpair))


                # . Inductive interactions by atom type
                # induct  alph   screen
                # 'H0'    0.0      1.0
                # (...)
                elif line.startswith ("induct"):
                    self.inductive = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, alpha, screen = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float, float])
                        inductive = EVBInductive (
                            evbType     = evbType   ,
                            alpha       = alpha     ,
                            screen      = screen    , )
                        self.inductive.append (inductive)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        ninductive = len (self.inductive)
                        npar += ninductive
                        print (messages["induct"] % (_MODULE_LABEL, ninductive))


                # . Inductive interactions by atom pair
                # a_induct   alph
                # 'C0' 'L-' 0.60
                # (...)
                elif line.startswith ("a_induct"):
                    self.ipairs = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        typea, typeb, alpha = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float])
                        ipair = EVBIndPair (
                            typea   =   typea   ,
                            typeb   =   typeb   ,
                            alpha   =   alpha   , )
                        self.ipairs.append (ipair)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nipair = len (self.ipairs)
                        npar += nipair
                        print (messages["a_induct"] % (_MODULE_LABEL, nipair))


                # . Electrostatic screening corrections by atom type
                # elect   mu_s
                # 'H0'    2.00
                # (...)
                elif line.startswith ("elect"):
                    self.screen = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        evbType, mus = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), float])
                        screen = EVBScreen (
                            evbType     =   evbType ,
                            mus         =   mus     , )
                        self.screen.append (screen)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nscreen = len (self.screen)
                        npar += nscreen
                        print (messages["elect"] % (_MODULE_LABEL, nscreen))


                # . Electrostatic screening corrections by atom pair
                # a_elect     mu_s
                # 'H0'  'H0'  4.00
                # (...)
                elif line.startswith ("a_elect"):
                    self.spairs = []
                    line, comment = self._GetLineWithComment (data)
                    while line != "":
                        typea, typeb, mus = TokenizeLine (line, converters=[lambda convert: convert.replace ("'", ""), lambda convert: convert.replace ("'", ""), float])
                        spair = EVBScreenPair (
                            typea       =   typea   ,
                            typeb       =   typeb   ,
                            mus         =   mus     , )
                        self.spairs.append (spair)
                        line, comment = self._GetLineWithComment (data)
                    if logging:
                        nspair = len (self.spairs)
                        npar += nspair
                        print (messages["a_elect"] % (_MODULE_LABEL, nspair))
        except StopIteration:
            pass
        # . Close the file
        data.close ()
        if logging:
            print ("# . %s> Total %d parameters read" % (_MODULE_LABEL, npar))


    def PurgeTypes (self, types):
        """Purge library leaving only selected atom types."""
        if hasattr (self, "pairs"):
            pairs = []
            for pair in self.pairs:
                if (pair.typea in types) and (pair.typeb in types):
                    pairs.append (pair)
            self.pairs = pairs
        if hasattr (self, "morse"):
            morse = []
            for bond in self.morse:
                if bond.evbType in types:
                    morse.append (bond)
            self.morse = morse
        if hasattr (self, "angles"):
            angles = []
            for angle in self.angles:
                if angle.evbType in types:
                    angles.append (angle)
            self.angles = angles
        if hasattr (self, "torsions"):
            torsions = []
            for torsion in self.torsions:
                if (torsion.typea in types) and (torsion.typeb in types):
                    torsions.append (torsion)
            self.torsions = torsions
        if hasattr (self, "impropers"):
            impropers = []
            for improper in self.impropers:
                if improper.evbType in types:
                    impropers.append (improper)
            self.impropers = impropers
        if hasattr (self, "exnonbond"):
            exnonbond = []
            for nonbond in self.exnonbond:
                if nonbond.evbType in types:
                    exnonbond.append (nonbond)
            self.exnonbond = exnonbond
        if hasattr (self, "xpairs"):
            xpairs = []
            for pair in self.xpairs:
                if (pair.typea in types) and (pair.typeb in types):
                    xpairs.append (pair)
            self.xpairs = xpairs
        if hasattr (self, "vdwevb"):
            vdwevb = []
            for vdw in self.vdwevb:
                if vdw.evbType in types:
                    vdwevb.append (vdw)
            self.vdwevb = vdwevb
        if hasattr (self, "vpairs"):
            vpairs = []
            for pair in self.vpairs:
                if (pair.typea in types) and (pair.typeb in types):
                    vpairs.append (pair)
            self.vpairs = vpairs
        if hasattr (self, "solvdw"):
            solvdw = []
            for vdw in self.solvdw:
                if vdw.evbType in types:
                    solvdw.append (vdw)
            self.solvdw = solvdw
        if hasattr (self, "inductive"):
            inductives = []
            for inductive in self.inductive:
                if inductive.evbType in types:
                    inductives.append (inductive)
            self.inductive = inductives
        if hasattr (self, "ipairs"):
            ipairs = []
            for ipair in self.ipairs:
                if (ipair.typea in types) and (ipair.typeb in types):
                    ipairs.append (ipair)
            self.ipairs = ipairs
        if hasattr (self, "screen"):
            screens = []
            for screen in self.screen:
                if screen.evbType in types:
                    screens.append (screen)
            self.screen = screens
        if hasattr (self, "spairs"):
            spairs = []
            for spair in self.spairs:
                if (spair.typea in types) and (spair.typeb in types):
                    spairs.append (spair)
            self.spairs = spairs


    def WriteLibrary (self, filename="", digits=1):
        """Write library to a file or stdout."""
        template = {
            "pairs"      :   "'%s'  '%s'   %X.Yf    %X.Yf    %X.Yf    %X.Yf    %X.Yf"  ,
            "morse"      :   "'%s'         %X.Yf    %X.Yf"  ,
            "angles"     :   "'%s'         %X.Yf    %X.Yf    %X.Yf    %X.Yf"  ,
            "torsions"   :   "'%s'  '%s'   %X.Yf    %X.Yf    %X.Yf"  ,
            "impropers"  :   "'%s'         %X.Yf    %X.Yf    %X.Yf"  ,
            "exnonbond"  :   "'%s'         %X.Yf    %X.Yf    %X.Yf    %X.Yf"  ,
            "xpairs"     :   "'%s'  '%s'   %X.Yf    %X.Yf    %X.Yf    %X.Yf"  ,
            "vdwevb"     :   "'%s'         %X.Yf    %X.Yf"  ,
            "vpairs"     :   "'%s'  '%s'   %X.Yf    %X.Yf"  ,
            "solvdw"     :   "'%s'         %X.Yf    %X.Yf"  ,
            "inductive"  :   "'%s'         %X.Yf    %X.Yf"  ,
            "ipairs"     :   "'%s'  '%s'   %X.Yf"   ,
            "screen"     :   "'%s'         %X.Yf"   ,
            "spairs"     :   "'%s'  '%s'   %X.Yf"   , }
        formats = {}
        for (key, string) in template.iteritems ():
            formats[key] = string.replace ("X", "%d" % (6 + digits)).replace ("Y", "%d" % digits)
        lines = []
        if hasattr (self, "pairs"):
            lines.append ("m_pair         morse_ab    r_ab      beta     f_harmonic   r0_harmonic")
            for pair in self.pairs:
                lines.append (formats["pairs"] % (pair.typea, pair.typeb, pair.morseAB, pair.rab, pair.beta, pair.forceHarmonic, pair.radiusHarmonic))
            lines.append ("")
        if hasattr (self, "morse"):
            lines.append ("morse_type     morse_d    radius")
            for bond in self.morse:
                lines.append (formats["morse"] % (bond.evbType, bond.morseD, bond.radius))
            lines.append ("")
        if hasattr (self, "angles"):
            lines.append ("angle_type     force      angle0")
            for angle in self.angles:
                lines.append (formats["angles"] % (angle.evbType, angle.force, angle.angle0, angle.foo, angle.bar))
            lines.append ("")
        if hasattr (self, "torsions"):
            lines.append ("torsion_type   force      periodicity  phase_angle")
            for torsion in self.torsions:
                lines.append (formats["torsions"] % (torsion.typea, torsion.typeb, torsion.force, torsion.periodicity, torsion.phase))
            lines.append ("")
        if hasattr (self, "impropers"):
            lines.append ("itorsion_type  force      periodicity  itorsion_angle0")
            for improper in self.impropers:
                lines.append (formats["impropers"] % (improper.evbType, improper.force, improper.periodicity, improper.angle0))
            lines.append ("")
        if hasattr (self, "exnonbond"):
            lines.append ("exnonbond_type  ex_repl   beta     q_vdwa    q_vdwb")
            for nonbond in self.exnonbond:
                lines.append (formats["exnonbond"] % (nonbond.evbType, nonbond.exRepulsive, nonbond.beta, nonbond.qvdwa, nonbond.qvdwb))
            lines.append ("")
        if hasattr (self, "xpairs"):
            lines.append ("x_pair          ex_repl    beta     q_vdwa    q_vdwb")
            for pair in self.xpairs:
                lines.append (formats["xpairs"] % (pair.typea, pair.typeb, pair.exRepulsive, pair.beta, pair.qvdwa, pair.qvdwb))
            lines.append ("")
        if hasattr (self, "vdwevb"):
            lines.append ("vdwevb              vdwa      vdwb")
            for vdw in self.vdwevb:
                lines.append (formats["vdwevb"] % (vdw.evbType, vdw.repulsive, vdw.attractive))
            lines.append ("")
        if hasattr (self, "vpairs"):
            lines.append ("v_pair              vdwa      vdwb")
            for pair in self.vpairs:
                lines.append (formats["vpairs"] % (pair.typea, pair.typeb, pair.repulsive, pair.attractive))
            lines.append ("")
        if hasattr (self, "solvdw"):
            lines.append ("solvdw              vdwa      vdwb")
            for vdw in self.solvdw:
                lines.append (formats["solvdw"] % (vdw.evbType, vdw.repulsive, vdw.attractive))
            lines.append ("")
        if hasattr (self, "inductive"):
            lines.append ("induct  alph   screen")
            for inductive in self.inductive:
                lines.append (formats["inductive"] % (inductive.evbType, inductive.alpha, inductive.screen))
            lines.append ("")
        if hasattr (self, "ipairs"):
            lines.append ("a_induct   alph")
            for ipair in self.ipairs:
                lines.append (formats["ipairs"] % (ipair.typea, ipair.typeb, ipair.alpha))
            lines.append ("")
        if hasattr (self, "screen"):
            lines.append ("elect   mu_s")
            for screen in self.screen:
                lines.append (formats["screen"] % (screen.evbType, screen.mus))
            lines.append ("")
        if hasattr (self, "spairs"):
            lines.append ("a_elect     mu_s")
            for spair in self.spairs:
                lines.append (formats["spairs"] % (spair.typea, spair.typeb, spair.mus))
            lines.append ("")
        if filename != "":
            fo = open (filename, "w")
            for line in lines:
                fo.write ("%s\n" % line)
            fo.close ()
        else:
            for line in lines:
                print line


    def GetSolVDW (self, atomType):
        if hasattr (self, "solvdw"):
            for vdw in self.solvdw:
                if vdw.evbType == atomType:
                    return (vdw.repulsive, vdw.attractive)
        return None


    def GetEvbVDW (self, atomType):
        if hasattr (self, "vdwevb"):
            for vdw in self.vdwevb:
                if vdw.evbType == atomType:
                    return (vdw.repulsive, vdw.attractive)
        return None


    def GetBond (self, typea, typeb):
        # . Parameters for Morse pairs take precedence
        if hasattr (self, "pairs"):
            for pair in self.pairs:
                if (pair.typea == typea and pair.typeb == typeb) or (pair.typea == typeb and pair.typeb == typea):
                    (morseD, r0) = pair.morseAB, pair.rab
                    return (morseD, r0)
        # . No pairs found, search one-atom parameters
        if hasattr (self, "morse"):
            morsea = None
            morseb = None
            for morse in self.morse:
                if   morse.evbType == typea:
                    morsea = morse
                elif morse.evbType == typeb:
                    morseb = morse
            if morsea and morseb:
                # . Apply combination rules
                morseD = math.sqrt (morsea.morseD * morseb.morseD)
                r0     = morsea.radius + morseb.radius
                return (morseD, r0)
        return None


    def GetAngle (self, typea, typeb, typec):
        if hasattr (self, "angles"):
            for angle in self.angles:
                if (angle.typeb == typeb):
                    return (angle.force, angle.angle0)
        return None


    def GetTorsion (self, typeb, typec):
        if hasattr (self, "torsions"):
            for torsion in self.torsions:
                if (torsion.typea == typeb and torsion.typeb == typec) or (torsion.typea == typec and torsion.typeb == typeb):
                    return (torsion.force, torsion.periodicity, torsion.phase)
        return None


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__":
    library  = EVBLibrary (logging=True)
    library.WriteLibrary ()
