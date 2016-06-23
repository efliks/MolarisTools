This example shows how to run QM(ai)/MM alchemical mutagenesis.

On going from state I to state II, two hydrogen atoms are replaced by dummy atoms 
and a carbon atom is replaced by an oxygen atom.

Thus, the two states differ in the number of atoms. These numbers have to be adjusted
in the call_gaussian_fep.py script by modifying lines starting from STATE_I and STATE_II.

Each of the two states will have its own set of Gaussian files (input, output, checkpoint).
The reason for separating these files is that we want to take advantage of restarting the 
calculation of orbitals from a previous QM run. Since the two states differ in the 
composition of atoms, use of the same checkpoint file would lead to a mess.

--------------
Mikolaj Feliks
Thu Jun 23 12:40:41 PDT 2016
