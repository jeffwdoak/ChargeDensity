#!/usr/bin/env python

from unitcell import *
import numpy as np

locname1 = str(sys.argv[1])
locname2 = str(sys.argv[2])
locfile1 = open(locname1,'r')
locfile2 = open(locname2,'r')
loc1 = UnitCell()
loc1.read_poscar(locfile1,vel=False)
loc1.read_chgcar(locfile1)
loc2 = UnitCell()
loc2.read_poscar(locfile2,vel=False)
loc2.read_chgcar(locfile2)
den_z1 = loc1.integrate_z_density()
den_z2 = loc2.integrate_z_density()
z_pos = np.linspace(0,loc1.cell_vec[2,2],len(den_z1))
spline1 = loc1.spline_density_z(den_z1)
spline2 = loc2.spline_density_z(den_z2)
z_diff,diff = loc1.dens_diff_z(spline1,spline2)
for i in range(len(diff)):
    print z_diff[i],diff[i]
