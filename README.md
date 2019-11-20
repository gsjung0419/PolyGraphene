# PolyGraphene
Generation of polycrystalline graphene models

A. Libraries should be prepared. The paths are described in the makefile
1. fftw 2.voro++ 3. lammps 

B. All code and libraries are only tested with intel compiler + openmpi in 2013. 
There might be version problems with the updated compilers
GNU compiler has not been tested

C. The required inputs are in the file “grains.cpp”’s constructor in the src directory. 
Lammps input is “relax.in"

D. There are some redundant files in the src. The main stuff is in “main.cpp” and “grain.cpp”

E. The example is for 6nm randomly distributed grains. 
rg6i.pdb is the initial geometry before refinement
rg6e.pdb is the final geometry after refinement

The code was not designed for the public use. (minimum comments without checking errors)

2019.11.20 Uploaded!
The code may be updated or a new version will be uploaded.

Please contact if you have questions gs4phone@gmail.com or jungg@ornl.gov
Cite: Jung, GangSeob, Zhao Qin, and Markus J. Buehler. "Molecular mechanics of polycrystalline graphene with enhanced fracture toughness." Extreme Mechanics Letters 2 (2015): 52-59. 
