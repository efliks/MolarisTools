#-------------------------------------------------------------------------------
# . File      : DetermineBAT.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from MolarisTools.Units    import DEFAULT_PARM_LIB, DEFAULT_AMINO_LIB
from MolarisTools.Parser   import MolarisOutputFile
from MolarisTools.Library  import AminoLibrary, ParametersLibrary

_DEFAULT_FORCE  = 5.


def DetermineBAT (fileLibrary=DEFAULT_AMINO_LIB, fileMolarisOutput="determine_atoms.out", residueLabels=(), fileParameters=DEFAULT_PARM_LIB):
    """Determine bonds, angles and torsions to analyze their statistical distributions."""

    # . Load the log file
    mof        = MolarisOutputFile (fileMolarisOutput)
    # . Load the library
    library    = AminoLibrary (fileLibrary, logging=False)
    # . Load Enzymix parameters
    parameters = ParametersLibrary (fileParameters) if fileParameters else None

    if mof.nresidues > 0:
        for residue in mof.residues:
            include = True
            if len (residueLabels) > 0:
                if residue.label not in residueLabels:
                    include = False
            if include:
                component = library[residue.label]
                component.GenerateAngles (logging=False)
                component.GenerateTorsions (logging=False)

                # . Write bonds
                bondTypes, bondUnique = component._BondsToTypes ()
                for (bonda, bondb), (typea, typeb) in zip (component.bonds, bondTypes):
                    serials = []
                    for bond in (bonda, bondb):
                        for atom in residue.atoms:
                            if atom.label == bond: break
                        serials.append (atom.serial)
                    if parameters:
                        parBond = parameters.GetBond (typea, typeb)
                        if parBond:
                            print ("constraint_pair  %4d  %4d    %4.1f    %5.2f    # %4s    %4s" % (serials[0], serials[1], _DEFAULT_FORCE, parBond.r0, bonda, bondb))
                        else:
                            print ("constraint_pair  %4d  %4d    %4.1f    XXXXX    # %4s    %4s" % (serials[0], serials[1], _DEFAULT_FORCE,             bonda, bondb))
                    else:
                        print ("%4d    %4d    # %4s    %4s" % (serials[0], serials[1], bonda, bondb))

                # . Write angles
                angleTypes, angleUnique = component._AnglesToTypes ()
                for (anglea, angleb, anglec), (typea, typeb, typec) in zip (component.angles, angleTypes):
                    serials = []
                    for angle in (anglea, angleb, anglec):
                        for atom in residue.atoms:
                            if atom.label == angle: break
                        serials.append (atom.serial)
                    if parameters:
                        parAngle = parameters.GetAngle (typea, typeb, typec)
                        if parAngle:
                            print ("constraint_ang   %4d  %4d  %4d    %4.1f    %5.2f    # %4s    %4s    %4s" % (serials[0], serials[1], serials[2], _DEFAULT_FORCE, parAngle.r0, anglea, angleb, anglec))
                        else:
                            print ("constraint_ang   %4d  %4d  %4d    %4.1f    XXXXX    # %4s    %4s    %4s" % (serials[0], serials[1], serials[2], _DEFAULT_FORCE,              anglea, angleb, anglec))
                    else:
                        print ("%4d    %4d    %4d    # %4s    %4s    %4s" % (serials[0], serials[1], serials[2], anglea, angleb, anglec))

                # . Write torsions
                torsionTypes, torsionUnique, torsionGeneral = component._TorsionsToTypes ()
                for (torsiona, torsionb, torsionc, torsiond), (typea, typeb, typec, typed) in zip (component.torsions, torsionTypes):
                    serials = []
                    for torsion in (torsiona, torsionb, torsionc, torsiond):
                        for atom in residue.atoms:
                            if atom.label == torsion: break
                        serials.append (atom.serial)
                    if parameters:
                        parTorsion = parameters.GetTorsion (typeb, typec)
                        if parTorsion:
                            pass
                    else:
                        print ("%4d    %4d    %4d    %4d    # %4s    %4s    %4s    %4s" % (serials[0], serials[1], serials[2], serials[3], torsiona, torsionb, torsionc, torsiond))


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
