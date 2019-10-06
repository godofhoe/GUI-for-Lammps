# GUI-for-Lammps
A GUI to integrate the energy of beams on a user-define curve

To readfile, using:
output_S1 = ReadAngle('1')
output_S = ReadAngle('2')
unrolled = ReadUnrolled('3')

use 1.2.txt to be '1' and '2', unrolled60150 to be '3'.

In these file, 
c_1 : bending 
c_2 : stretching

In unrolled60150, the lines below "Atoms" indicate:
id 0 1 x y z
