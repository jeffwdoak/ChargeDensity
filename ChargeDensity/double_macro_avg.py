#!/usr/bin/env python

# double_macro_avg.py v1.0 9-19-2012 Jeff Doak jeff.w.doak@gmail.com

from chargedensity import *
import numpy as np
import sys

p1 = float(sys.argv[1])
p2 = float(sys.argv[2])
if len(sys.argv) > 3:
    if str(sys.argv[3]) == "CHG" :
        a = ChargeDensity(str(sys.argv[3]),format_="chgcar")
    else:
        a = ChargeDensity(str(sys.argv[3]))
else:
    a = ChargeDensity("LOCPOT")
den_z = a.integrate_z_density()
z_pos = np.linspace(0,a.unitcell.cell_vec[2,2],len(den_z))
macro_z = a.macro_mismatch(p1,p2)
for i in range(len(den_z)):
    print z_pos[i],den_z[i],macro_z[i]
