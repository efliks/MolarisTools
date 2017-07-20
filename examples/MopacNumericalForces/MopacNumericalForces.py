#-------------------------------------------------------------------------------
# . File      : MopacNumericalForces.py
# . Program   : MolarisTools
# . Copyright : USC, Mikolaj Feliks (2015-2017)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, copy, collections, subprocess

from MolarisTools.Parser  import MopacInputFile, MopacOutputFile


Atom = collections.namedtuple ("Atom" , "label x y z")

mopac      =  os.path.join (os.environ["HOME"], "local", "bin", "MOPAC2009.exe")
inputFile  =  MopacInputFile ()
delta      =  0.05
workdir    =  "work"

# . Modify keywords
# inputFile.keywords.remove ("MOL_QMMM")


print ("\n=== Numerical derivatives ===")

for (atomSerial, atom) in enumerate (inputFile.atoms, 1):
    derivatives = (
        (
            (
            atom.x  +  delta  ,
            atom.y            ,
            atom.z            ,
            os.path.join (workdir, "atom%02d_xp.mop" % atomSerial) ,
            ) ,

            (
            atom.x  -  delta  ,
            atom.y            ,
            atom.z            ,
            os.path.join (workdir, "atom%02d_xm.mop" % atomSerial) ,
            ) ,
        ),

        (
            (
            atom.x            ,
            atom.y  +  delta  ,
            atom.z            ,
            os.path.join (workdir, "atom%02d_yp.mop" % atomSerial) ,
            ) ,

            (
            atom.x            ,
            atom.y  -  delta  ,
            atom.z            ,
            os.path.join (workdir, "atom%02d_ym.mop" % atomSerial) ,
            ) ,
        ),

        (
            (
            atom.x            ,
            atom.y            ,
            atom.z  +  delta  ,
            os.path.join (workdir, "atom%02d_zp.mop" % atomSerial) ,
            ) ,

            (
            atom.x            ,
            atom.y            ,
            atom.z  -  delta  ,
            os.path.join (workdir, "atom%02d_zm.mop" % atomSerial) ,
            ) ,
        ),
    )

    results = []
    for derivative in derivatives:
        # . Run calculations
        energies = []
        for direction in derivative:
            (x, y, z, fileName) = direction

            newFile   = copy.deepcopy (inputFile)
            atomIndex = atomSerial - 1
            newAtom = Atom (
                label   =   atom.label  ,
                x       =   x           ,
                y       =   y           ,
                z       =   z           ,
                )
            newFile.atoms[atomIndex] = newAtom
            newFile.Write (fileName)
            # print ("Wrote file %s." % fileName)

            # . Run MOPAC
            stem, extension = os.path.splitext (fileName)
            logFile = stem + ".out"
            if not os.path.exists (logFile):
                subprocess.check_call ([mopac, fileName])

            # . Read output file
            mof = MopacOutputFile (logFile)
            energies.append (mof.Efinal) 


        # . Analyze results
        (forward, reverse) = derivative
        (forEnergy, revEnergy) = energies

        derv = -1. * (forEnergy - revEnergy) / (2. * delta)
        results.append (derv)


    # . For each atom, print its numerical derivatives
    (fx, fy, fz) = results
    print ("%2s     %9.4f     %9.4f     %9.4f" % (atom.label, fx, fy, fz))



print ("\n=== Analytic derivatives ===")
fileAnalytic = MopacOutputFile ("run.out")

for (atom, force) in zip (inputFile.atoms, fileAnalytic.forces):
    print ("%2s     %9.4f     %9.4f     %9.4f" % (atom.label, force.x, force.y, force.z))
