#!/usr/bin/env python

# plot_chg_den_z.py v1.0 9-20-2012 Jeff Doak jeff.w.doak@gmail.com

from chargedensity import *
import numpy as np

a = ChargeDensity("CHG",format_="chgcar")
den_z = a.integrate_z_density()
z_pos = np.linspace(0,a.unitcell.cell_vec[2,2],len(den_z))
print len(den_z)
for i in range(len(den_z)):
    print z_pos[i],den_z[i]
