The following example demonstrates how to run a simulation in
Molaris with a QM(ai)/MM potential on two states.

Although the potential used here is the same for both states,
any combination of QM(ai)/MM potentials can be used (however, the Python
script would have to be slightly modified).

The difference between the two states is that atoms 9 and 10 are
alchemically transformed from hydrogens to fluorines.

The Python script runs calculations for both states simultaneously,
which is a 2X improvement in speed compared to how it was done 
previously in Molaris (sequentially).

Note that the simulation now requires twice as many CPUs (16 instead of
the usual 8).

To run the example:
    (1) Generate input files: python prep_equil.py

    (2) Allocate the job (on the San Diego cluster, comet): sbatch < alloc_regular_node.bash


This example was tested with molaris_hpc9.15, whose md5sum was 
61a84fb705d3c71dba8267ecd2990ff0.


--------------
Mikolaj Feliks
Fri Jul 22 11:22:24 PDT 2016
