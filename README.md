# MolarisTools
A Python toolkit to facilitate working with Molaris-XG.

Author: Mikolaj Feliks <<mikolaj.feliks@gmail.com>><br>
Released as open source software under the GNU GPL v3.0 license (see COPYING).


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
be obtained from the website of the Warshel Group:
<<http://laetro.usc.edu/software.html>>

To install MolarisTools, first download and unpack the ZIP archive
from this website. Prepare the package for installation:

python setup.py sdist

In the next step, build and install the package:

python setup.py build<br>
python setup.py install --prefix=${HOME}/local

The "prefix" argument is optional and indicates where MolarisTools 
should be placed. I keep a "local" directory in my home 
directory where all Python packages end up.


MolarisTools are actively developed and new features are added 
as they are needed. Major changes to the code can happen anytime.
Bugs are inevitable - please report if you have found one.
