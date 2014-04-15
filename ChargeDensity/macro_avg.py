#!/usr/bin/env python

# macro_avg.py v1.0 9-19-2012 Jeff Doak jeff.w.doak@gmail.com

from chargedensity import *
import numpy as np
import sys

if len(sys.argv) > 1:
    if str(sys.argv[1]) == "CHG":
        a = ChargeDensity(str(sys.argv[1]),format_="chgcar")
    else:
        a = ChargeDensity(str(sys.argv[1]))
else:
    a = ChargeDensity("LOCPOT")
if len(sys.argv) > 2:
    p1 = float(sys.argv[2])
else:
    p1 = a.unitcell.cell_vec[2,2]/float(10)
den_z = a.integrate_z_density()
z_pos = np.linspace(0,a.unitcell.cell_vec[2,2],len(den_z))
macro_z = a.macro_avg_z(p1)
for i in range(len(den_z)):
    print z_pos[i],den_z[i],macro_z[i]
