The following example shows how to run QM/MM molecular dynamics 
in Molaris with increasing temperature.


To run the simulation:
----------------------

(1) Edit heating_template.inp to choose Mopac (semi-empirical QM/MM) 
or Gaussian (ab initio QM/MM). Comment out/in relevant lines.

With Gaussian, you will have a choice to run either mechanical 
embedding (ME) or electrostatic embedding (EE).

ME is a crude approach to QM/MM because it does not allow for 
the polarization of the QM region.

In ME, electrostatic interactions between the MM and QM regions are 
calculated at the classical level.

The default is to use EE, where the QM region can be polarized. See the 
comments inside heating_template.inp and CallGaussian.py on how to 
run EE or ME.


(2) Generate input files for heating: python PrepareInputFiles.py


(3) Run the simulation: bash run_all.bash


Keep in mind that simulations with Gaussian are very slow compared
to Mopac or EVB if only one CPU is used.

Molaris libraries used in this example are located in libs/. Note that 
pdb_dictionary and solvent.opt are not included. These should be
taken from the official distribution of Molaris.


This example was tested with molaris_hpc9.15, whose md5sum
was 612359cd225fa825600203ede8da76bb (because version numbers are
meaningless in Molaris).
