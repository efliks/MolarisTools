#-------------------------------------------------------------------------------
# . File      : GenerateComponent.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import copy

from MolarisTools.Library  import AminoLibrary, AminoAtom


library = AminoLibrary (logging=False, filename="template.lib")

analogues = (
("EMH" ,  "Analogue -CF2- " ,   "EMF" ,  (("C3B" , "C3B" , "CT" ,  0.50) ,   ("H1" , "F1"  , "F-" , -0.20) ,   ("H2" , "F2"  , "F-" , -0.20))  ) ,
("EMH" ,  "Analogue -CCl2-" ,   "EMC" ,  (("C3B" , "C3B" , "CT" ,  0.50) ,   ("H1" , "CL1" , "CL" , -0.20) ,   ("H2" , "CL2" , "CL" , -0.20))  ) ,
("EMH" ,  "Analogue -CBr2-" ,   "EMB" ,  (("C3B" , "C3B" , "CT" ,  0.50) ,   ("H1" , "BR1" , "BR" , -0.20) ,   ("H2" , "BR2" , "BR" , -0.20))  ) ,
)

for ianalogue, (templateLabel, title, label, atoms) in enumerate (analogues, 1):
    component    = library[templateLabel]

    # . Clone residue
    clone        = copy.deepcopy (component)
    clone.serial = ianalogue + 100
    clone.name   = label

    # . Replace atoms
    for (atomLabel, atomLabelNew, atomTypeNew, atomChargeNew) in atoms:
        clone.ReplaceAtom (atomLabel, atomLabelNew, atomTypeNew, atomChargeNew)

    # . Write new residue
    clone.Write (title="%s" % title, showGroups=True, showLabels=True, terminate=True, filename="analogue_%s.lib" % label)
