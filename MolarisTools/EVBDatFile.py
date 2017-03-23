#-------------------------------------------------------------------------------
# . File      : EVBDatFile.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from    Utilities         import TokenizeLine
from    Units             import DEFAULT_EVB_LIB
from    MolarisInputFile  import MolarisInputFile
import  collections, exceptions, os, math

EVBContainer = collections.namedtuple ("EVBContainer", "serials  types  exist")
EVBBond     = collections.namedtuple ("EVBBond"    , "r0  alpha  diss  f_harm  r_harm  ptype")
EVBAngle    = collections.namedtuple ("EVBAngle"   , "angle0  force  gausD  gausig  ptype")
EVBTorsion  = collections.namedtuple ("EVBTorsion" , "force  n  phase  ptype")


_DEFAULT_FILENAME = os.path.join ("evb_heat_01", "evb.dat")
_MODULE_LABEL     = "EVBDat"
_ATOMS_PER_LINE   = 10
_VMD_DIR          = "vmd"


class EVBDatFile (object):
    """A class to represent an evb.dat file."""

    def __init__ (self, filename=_DEFAULT_FILENAME, logging=True):
        """Constructor."""
        self.filename = filename
        self._Parse (logging=logging)


    def _GetContainer (self, data, ntokens):
        serials   = TokenizeLine (next (data), converters=([int, ] * ntokens    ))
        types     = TokenizeLine (next (data), converters=([int, ] * self.nforms))
        exist     = TokenizeLine (next (data), converters=([int, ] * self.nforms))
        container = EVBContainer (
            serials =   serials ,
            types   =   types   ,
            exist   =   exist   , )
        return container


    def _Parse (self, logging):
        data = open (self.filename, "r")
        if logging:
            print ("# . %s> Parsing file \"%s\"" % (_MODULE_LABEL, self.filename))
        try:
            while True:
                line = next (data)
                if line.count ("# of evb atoms, # of resforms"):
                    (natoms, nforms) = TokenizeLine (line, converters=[int, int])
                    nlines  = (natoms / _ATOMS_PER_LINE) + (1 if (natoms % _ATOMS_PER_LINE) > 0 else 0)
                    # . Read atom serial numbers
                    serials = []
                    for i in range (nlines):
                        tokens = TokenizeLine (next (data), converters=([int, ] * _ATOMS_PER_LINE))
                        for token in tokens:
                            if token != None:
                                serials.append (token)
                    # . Read atom type numbers for each form
                    forms = []
                    for i in range (nforms):
                        types = []
                        for j in range (nlines):
                            tokens = TokenizeLine (next (data), converters=([int, ] * _ATOMS_PER_LINE))
                            for token in tokens:
                                if token != None:
                                    types.append (token)
                        forms.append (types)
                    self.types   = forms
                    # . Read atomic charges for each form
                    charges = []
                    for i in range (nforms):
                        form = []
                        for j in range (nlines):
                            tokens = TokenizeLine (next (data), converters=([float, ] * _ATOMS_PER_LINE))
                            for token in tokens:
                                if token != None:
                                    form.append (token)
                        charges.append (form)
                    self.charges = charges
                    self.natoms  = natoms
                    self.nforms  = nforms
                    if logging:
                        print ("# . %s> Found %d EVB atoms in %d forms" % (_MODULE_LABEL, natoms, nforms))
                        types = set ()
                        for form in self.types:
                            for atom in form:
                                types.add (atom)
                        ntypes = len (types)
                        print ("# . %s> Found %d unique atom types" % (_MODULE_LABEL, ntypes))

                elif line.count ("bonds(atoms,types,exist)"):
                    (nbonds, foo) = TokenizeLine (line, converters=[int, None])
                    bonds = []
                    for i in range (nbonds):
                        bond = self._GetContainer (data, 2)
                        bonds.append (bond)
                    self.bonds = bonds
                    if logging:
                        print ("# . %s> Found %d EVB bonds" % (_MODULE_LABEL, nbonds))

                elif line.count ("angles(atoms,types,exist)"):
                    (nangles, foo) = TokenizeLine (line, converters=[int, None])
                    angles = []
                    for i in range (nangles):
                        angle = self._GetContainer (data, 3)
                        angles.append (angle)
                    self.angles = angles
                    if logging:
                        print ("# . %s> Found %d EVB angles" % (_MODULE_LABEL, nangles))

                elif line.count (" torsions(atoms,types,exist)"):
                    (ntorsions, foo) = TokenizeLine (line, converters=[int, None])
                    torsions = []
                    for i in range (ntorsions):
                        torsion = self._GetContainer (data, 4)
                        torsions.append (torsion)
                    self.torsions = torsions
                    if logging:
                        print ("# . %s> Found %d EVB torsions" % (_MODULE_LABEL, ntorsions))

                elif line.count ("morse potential parameters(r0,alpha,diss,f_harm,r_harm,type)"):
                    (nparameters, foo) = TokenizeLine (line, converters=[int, None])
                    self.parBonds = []
                    for i in range (nparameters):
                        (r0, alpha, diss, f_harm, r_harm, ptype) = TokenizeLine (next (data), converters=[float, float, float, float, float, int])
                        bond = EVBBond (
                            r0      =   r0      ,
                            diss    =   diss    ,
                            alpha   =   alpha   ,
                            f_harm  =   f_harm  ,
                            r_harm  =   r_harm  ,
                            ptype   =   ptype   ,  )
                        self.parBonds.append (bond)
                    if logging:
                        print ("# . %s> Found %d EVB Morse parameters" % (_MODULE_LABEL, nparameters))

                elif line.count ("angle parameters(angle0(radian),force,gausD,gausig,type"):
                    (nparameters, foo) = TokenizeLine (line, converters=[int, None])
                    self.parAngles = []
                    for i in range (nparameters):
                        (angle0, force, gausD, gausig, ptype) = TokenizeLine (next (data), converters=[float, float, float, float, int])
                        angle = EVBAngle (
                            angle0  =   angle0  ,
                            force   =   force   ,
                            gausD   =   gausD   ,
                            gausig  =   gausig  ,
                            ptype   =   ptype   , )
                        self.parAngles.append (angle)
                    if logging:
                        print ("# . %s> Found %d EVB angle parameters" % (_MODULE_LABEL, nparameters))

                elif line.count ("torsion parameters(force,n,phase_angle,type)"):
                    if not line.count ("itorsion"):
                        (nparameters, foo) = TokenizeLine (line, converters=[int, None])
                        self.parTorsions = []
                        for i in range (nparameters):
                            (force, n, phase, ptype) = TokenizeLine (next (data), converters=[float, float, float, int])
                            torsion = EVBTorsion (
                                force   =   force   ,
                                n       =   n       ,
                                phase   =   phase   ,
                                ptype   =   ptype   , )
                            self.parTorsions.append (torsion)
                        if logging:
                            print ("# . %s> Found %d EVB torsion parameters" % (_MODULE_LABEL, nparameters))
        except StopIteration:
            pass
        # . File closing
        data.close ()


    def _GetConvert (self, filenameInput, state, logging, purify=False):
        convert = {}
        if (filenameInput != ""):
            mif = MolarisInputFile (filenameInput, logging=logging)
            for (i, atom) in enumerate (mif.states[(state - 1)], 1):
                (enzymixCharge, enzymixType, label, group) = atom.SplitComment ()
                if purify:
                    label = label.replace ("'", "p").replace ("*", "p")
                convert[i] = (label, atom.atype)
        return convert


    def Decode (self, filenameInput, state=1, digits=1, showOnly=(), extended=False, logging=True):
        """Translate a DAT file into a list of parameters."""
        template = {
            "bonds"         :   "%4s  %4s  %X.Yf    %X.Yf    %X.Yf    %X.Yf    %X.Yf      (%d)"  ,
            "angles"        :   "%4s  %4s  %4s  %X.Yf    %X.Yf    %X.Yf    %X.Yf      (%d)"  ,
            "torsions"      :   "%4s  %4s  %4s  %4s   %X.Yf    %X.Yf    %X.Yf      (%d)"  ,
            "bonds_ext"     :   "%4s  %4s  %X.Yf    %X.Yf    %X.Yf    %X.Yf    %X.Yf      (%d)   %s    %s"  ,
            "angles_ext"    :   "%4s  %4s  %4s  %X.Yf    %X.Yf    %X.Yf    %X.Yf      (%d)   %s    %s    %s"  ,
            "torsions_ext"  :   "%4s  %4s  %4s  %4s   %X.Yf    %X.Yf    %X.Yf      (%d)   %s    %s    %s    %s"  , }
        formats = {}
        for (key, string) in template.iteritems ():
            formats[key] = string.replace ("X", "%d" % (6 + digits)).replace ("Y", "%d" % digits)
        convert = self._GetConvert (filenameInput, state, logging, purify=False)
    
        print ("\n-- EVB bonds --")
        for bond in self.bonds:
            if bond.exist[(state - 1)]:
                found = True
                try:
                    (i, j) = bond.serials
                    (labela, typea), (labelb, typeb) = convert[i], convert[j]
                except:
                    found = False
                    if logging:
                        (seriala, serialb) = bond.serials
                        print ("# . %s> Warning: Bond (%d, %d) involves non-EVB atoms" % (_MODULE_LABEL, seriala, serialb))
                if found:
                    select   = bond.types[(state - 1)]
                    parBond  = self.parBonds[(select - 1)]
                    if showOnly:
                        if parBond.ptype not in showOnly:
                            continue
                    if extended:
                        print (formats["bonds_ext"] % (labela, labelb, parBond.r0, parBond.alpha, parBond.diss, parBond.f_harm, parBond.r_harm, parBond.ptype, typea, typeb))
                    else:
                        print (formats["bonds"]     % (labela, labelb, parBond.r0, parBond.alpha, parBond.diss, parBond.f_harm, parBond.r_harm, parBond.ptype))
    
        print ("\n-- EVB angles --")
        for angle in self.angles:
            if angle.exist[(state - 1)]:
                found = True
                try:
                    (i, j, k) = angle.serials
                    (labela, typea), (labelb, typeb), (labelc, typec) = convert[i], convert[j], convert[k]
                except:
                    found = False
                    if logging:
                        (seriala, serialb, serialc) = angle.serials
                        print ("# . %s> Warning: Angle (%d, %d, %d) involves non-EVB atoms" % (_MODULE_LABEL, seriala, serialb, serialc))
                if found:
                    select   = angle.types[(state - 1)]
                    parAngle = self.parAngles[(select - 1)]
                    if showOnly:
                        if parAngle.ptype not in showOnly:
                            continue
                    if extended:
                        print (formats["angles_ext"] % (labela, labelb, labelc, (parAngle.angle0 * 180. / math.pi), parAngle.force, parAngle.gausD, parAngle.gausig, parAngle.ptype, typea, typeb, typec))
                    else:
                        print (formats["angles"]     % (labela, labelb, labelc, (parAngle.angle0 * 180. / math.pi), parAngle.force, parAngle.gausD, parAngle.gausig, parAngle.ptype))
    
        print ("\n-- EVB torsions --")
        for torsion in self.torsions:
            if torsion.exist[(state - 1)]:
                found = True
                try:
                    (i, j, k, l) = torsion.serials
                    (labela, typea), (labelb, typeb), (labelc, typec), (labeld, typed) = convert[i], convert[j], convert[k], convert[l]
                except:
                    found = False
                    if logging:
                        (seriala, serialb, serialc, seriald) = torsion.serials
                        print ("# . %s> Warning: Torsion (%d, %d, %d, %d) involves non-EVB atoms" % (_MODULE_LABEL, seriala, serialb, serialc, seriald))
                if found:
                    select     = torsion.types[(state - 1)]
                    parTorsion = self.parTorsions[(select - 1)]
                    if showOnly:
                        if parTorsion.ptype not in showOnly:
                            continue
                    if extended:
                        print (formats["torsions_ext"] % (labela, labelb, labelc, labeld, parTorsion.force, parTorsion.n, (parTorsion.phase * 180. / math.pi), parTorsion.ptype, typea, typeb, typec, typed))
                    else:
                        print (formats["torsions"]     % (labela, labelb, labelc, labeld, parTorsion.force, parTorsion.n, (parTorsion.phase * 180. / math.pi), parTorsion.ptype))


    def GenerateVMDCommands (self, location=_VMD_DIR, state=1, filenameInput="", useLabels=False, logging=True):
        """Generate a set of VMD commands for measuring bonds, angles and torsions in a QM/MM trajectory."""
        convert  = self._GetConvert (filenameInput, state, logging, purify=useLabels)
        commands = []

        for (counter, bond) in enumerate (self.bonds):
            if bond.exist[(state - 1)]:
                if convert != {}:
                    try:
                        (i, j) = bond.serials
                        (labela, typea), (labelb, typeb) = convert[i], convert[j]
                        commands.append ("#  %s--%s" % (labela, labelb))
                    except:
                        commands.append ("# Bond (%d, %d) involves non-EVB atoms" % (seriala, serialb))
                (indexa, indexb) = map (lambda q: (q - 1), bond.serials)
                commands.append ("label  add    Bonds       0/%d     0/%d"        % (indexa, indexb))
                (seriala, serialb) = bond.serials
                if useLabels:
                    commands.append ("label  graph  Bonds       %3d     %s/dist_%s_%s.dat" % (counter, location, labela, labelb))
                else:
                    commands.append ("label  graph  Bonds       %3d     %s/dist_%d_%d.dat" % (counter, location, seriala, serialb))

        for (counter, angle) in enumerate (self.angles):
            if angle.exist[(state -1)]:
                if convert != {}:
                    try:
                        (i, j, k) = angle.serials
                        (labela, typea), (labelb, typeb), (labelc, typec) = convert[i], convert[j], convert[k]
                        commands.append ("#  %s--%s--%s" % (labela, labelb, labelc))
                    except:
                        commands.append ("# Angle (%d, %d, %d) involves non-EVB atoms" % (seriala, serialb, serialc))
                (indexa, indexb, indexc) = map (lambda q: (q - 1), angle.serials)
                commands.append ("label  add    Angles      0/%d     0/%d   0/%d"   % (indexa, indexb, indexc))
                (seriala, serialb, serialc) = angle.serials
                if useLabels:
                    commands.append ("label  graph  Angles      %3d     %s/angl_%s_%s_%s.dat" % (counter, location, labela, labelb, labelc))
                else:
                    commands.append ("label  graph  Angles      %3d     %s/angl_%d_%d_%d.dat" % (counter, location, seriala, serialb, serialc))

        for (counter, torsion) in enumerate (self.torsions):
            if torsion.exist[(state - 1)]:
                if convert != {}:
                    try:
                        (i, j, k, l) = torsion.serials
                        (labela, typea), (labelb, typeb), (labelc, typec), (labeld, typed) = convert[i], convert[j], convert[k], convert[l]
                        commands.append ("#  %s--%s--%s--%s" % (labela, labelb, labelc, labeld))
                    except:
                        commands.append ("# Torsion (%d, %d, %d, %d) involves non-EVB atoms" % (seriala, serialb, serialc, seriald))
                (indexa, indexb, indexc, indexd) = map (lambda q: (q - 1), torsion.serials)
                commands.append ("label  add    Dihedrals   0/%d     0/%d   0/%d   0/%d" % (indexa, indexb, indexc, indexd))
                (seriala, serialb, serialc, seriald) = torsion.serials
                if useLabels:
                    commands.append ("label  graph  Dihedrals   %3d     %s/dihe_%s_%s_%s_%s.dat"    % (counter, location, labela, labelb, labelc, labeld))
                else:
                    commands.append ("label  graph  Dihedrals   %3d     %s/dihe_%d_%d_%d_%d.dat"    % (counter, location, seriala, serialb, serialc, seriald))
        return commands


    def _ParseVMDFile (self, filename, wrapTorsions=False):
        lines   = open (filename, "r").readlines ()
        collect = []
        for line in lines:
            (foo, value) = TokenizeLine (line, converters=[float, float])
            if wrapTorsions:
                if (value < 0.):
                    if (value < 90.):
                        value += 360.
            collect.append (value)
        return collect


    def ReadVMDFiles (self, location=_VMD_DIR, logging=True, wrapTorsions=False):
        """Read .dat files generated by VMD.

        For now, works only with files described by serials, not labels."""
        if hasattr (self, "bonds"):
            self.bondsValues   = []
            for bond in self.bonds:
                (seriala, serialb) = bond.serials
                filename  = os.path.join (location, "dist_%d_%d.dat" % (seriala, serialb))
                if os.path.exists (filename):
                    if logging:
                        print ("# . %s> Parsing file %s..." % (_MODULE_LABEL, filename))
                    values  = self._ParseVMDFile (filename)
                    self.bondsValues.append (values)

        if hasattr (self, "angles"):
            self.anglesValues   = []
            for angle in self.angles:
                (seriala, serialb, serialc) = angle.serials
                filename  = os.path.join (location, "angl_%d_%d_%d.dat" % (seriala, serialb, serialc))
                if os.path.exists (filename):
                    if logging:
                        print ("# . %s> Parsing file %s..." % (_MODULE_LABEL, filename))
                    values  = self._ParseVMDFile (filename)
                    self.anglesValues.append (values)

        if hasattr (self, "torsions"):
            self.torsionsValues = []
            for torsion in self.torsions:
                (seriala, serialb, serialc, seriald) = torsion.serials
                filename  = os.path.join (location, "dihe_%d_%d_%d_%d.dat" % (seriala, serialb, serialc, seriald))
                if os.path.exists (filename):
                    if logging:
                        print ("# . %s> Parsing file %s..." % (_MODULE_LABEL, filename))
                    values  = self._ParseVMDFile (filename, wrapTorsions=wrapTorsions)
                    self.torsionsValues.append (values)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
