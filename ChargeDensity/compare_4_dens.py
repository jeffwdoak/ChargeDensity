#!/usr/bin/env python

# compare_4_dens.py v1.0 9-27-2012 Jeff Doak jeff.w.doak@gmail.com

from chargedensity import *
import numpy as np

# The first file should have the density grid that will be used for the final
# density that is returned.
locname1 = str(sys.argv[1])
locname2 = str(sys.argv[2])
locname3 = str(sys.argv[3])
locname4 = str(sys.argv[4])

locfile1 = open(locname1,'r')
locfile2 = open(locname2,'r')
locfile3 = open(locname3,'r')
locfile4 = open(locname4,'r')

loc1 = ChargeDensity(locfile1,"chgcar")
loc2 = ChargeDensity(locfile2,"chgcar")
loc3 = ChargeDensity(locfile3,"chgcar")
loc4 = ChargeDensity(locfile4,"chgcar")

den_z1 = loc1.integrate_z_density()
den_z2 = loc2.integrate_z_density()
den_z3 = loc3.integrate_z_density()
den_z4 = loc4.integrate_z_density()


z_pos = np.linspace(0,loc1.unitcell.cell_vec[2,2],len(den_z1))
spline1 = loc1.spline_density_z(den_z1)
spline2 = loc2.spline_density_z(den_z2)
spline3 = loc2.spline_density_z(den_z3)
spline4 = loc2.spline_density_z(den_z4)
z_diff,diff1 = loc1.dens_diff_z(spline1,spline2)
z_diff2,diff2 = loc1.dens_diff_z(spline3,spline4)
for i in range(len(diff1)):
    print z_diff[i],(diff1[i]-diff2[i])
