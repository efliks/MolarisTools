:warning: **IMPORTANT NOTE** :warning:

**This is a legacy library that is only compatible with Python 2.
Is is kept here solely for historical reasons.**


# MolarisTools
A Python toolkit to facilitate working with Molaris-XG.


_Key features:_
  * QM/MM interface with electrostatic embedding to Gaussian, Mopac, GAMESS-US, ORCA, Q-Chem
  * Parsing of Molaris files (input, log, gap, FVX, mol.in, evb.dat)
  * Parsing of files from quantum chemical packages (Gaussian, Mopac, GAMESS-US, ORCA, Q-Chem)
  * Parsing of geometry files (PDB, xyz, xyz trajectories)
  * Reading and writing of Molaris libraries (amino-library, ENZYMIX \& EVB parameters)
  * Handling of amino-components (calculation of partial charges, generation of angles and dihedrals, topology operations, merging)
  * Conversion between Molaris and CHARMM topology formats
  * Generation of tables for input files with EVB atoms \& bonds
  * Automatic generation of amino-components from PDB files based on coordinates and distances
  * Parsing of 1D \& 2D PES scans
  * LRA calculations


_Installation instructions:_

MolarisTools is a stand-alone Python library and as such does not 
need Molaris to be preinstalled. Nevertheless, a copy of Molaris can
be obtained from the [Warshel Group](http://laetro.usc.edu/software.html).

To install MolarisTools, clone the repository from GitHub (assuming that you have 
git installed on your computer):

```
git clone https://github.com/mfx9/MolarisTools.git
```

Or download and unpack the ZIP package from this website. In the next
step, adjust the PYTHONPATH variable so it points to the location
of MolarisTools, for example (in Bash):

```
export PYTHONPATH=${HOME}/MolarisTools:${PYTHONPATH}
```

Add the above line to your ~/.profile or ~/.bashrc file.

MolarisTools are actively developed and new features are added 
as they are needed. Major changes to the code can happen anytime.
Bugs are inevitable - please report if you have found one.
