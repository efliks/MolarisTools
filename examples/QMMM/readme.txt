To run QM/MM molecular dynamics in Molaris:

(1) Edit heating_template.inp to choose Mopac (semi-empirical QM/MM) or Gaussian (ab initio QM/MM).

With Gaussian, you will have a choice to run either mechanical embedding (ME) or 
electrostatic embedding (EE).

ME is a crude approach to QM/MM because it does not allow for the polarization of the QM region.
In ME, electrostatic interactions between the MM and QM regions are calculated at the classical level.

The default is to use EE, where the QM region can be polarized. See the comments inside heating_template.inp 
and call_gaussian.py on how to run EE or ME.


(2) Generate input files for heating: python prep_heating2.py


(3) Run the simulation: bash run_all.bash


This example was tested with molaris_hpc9.15, whose md5sum was 
0cf32647d1aa4087689b14a2c6d2220e (because version numbers are not terribly meaningful in Molaris).

Note that before running this example, you need to supply your own
libraries for parm.lib, pdb_dictionary and solvent.opt. 
They are part of the Molaris package and hence cannot be attached here.

--------------
Mikolaj Feliks
Fri May 20 15:37:19 PDT 2016
