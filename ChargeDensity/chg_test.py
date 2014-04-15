#!/usr/bin/env python

from unitcell import *
import numpy as np

loc = open("LOCPOT",'r')
a = UnitCell()
a.read_poscar(loc,vel=False)
a.read_chgcar(loc)
den_z = a.integrate_z_density()
z_pos = np.linspace(0,a.cell_vec[2,2],len(den_z))
lat_const = 4./float(len(a.generate_z_list()))
#lat_const = 10./a.cell_vec[-1,-1]
#print lat_const
#macro_z,macro_z2 = a.macro_avg_density_z(lat_const)
p1 = 6.557683331
p2 = 5.85347097
macro_z = a.macro_mismatch(p1,p2)
#macro_z = a.macro_avg_density_z()
print len(den_z)
for i in range(len(den_z)):
    print z_pos[i],den_z[i]
