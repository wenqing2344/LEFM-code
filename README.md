This is a demo of fracture in alloy system Nb45Ta25Ti15Hf15, a single crystal with fracture surface 112.
Structral data for publication is too big to be included on Github. 

1: To run this fracture simulation code, LAMMPS software with ACE machine learning potential library is required. 
The details of compiling LAMMPS with ACE can be found in this link: 
https://pacemaker.readthedocs.io/en/latest/pacemaker/quickstart/#lammps

2: Before running these codes, a few folders are needed. Please do the following
mkdir equi
mkdir boundary_displace
mkdir inner_displace
mkdir other_equi
mkdir dump

3: Open get_displacement.ipynb, run it to get the initial LEFM loaded simulation cell.

4: Run loop_relax.in with LAMMPS. The interatomic potential is f940.yace 

5: Fracture configurations are in folder equi. We expect a twin formation at the crack tip.
