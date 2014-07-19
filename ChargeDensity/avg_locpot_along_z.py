#!/usr/bin/python

# avg_locpot_along_z.py v1.0 9-18-2012 Jeff Doak jeff.w.doak@gmail.com

from chargedensity import *
import numpy as np
import sys

locfile = str(sys.argv[1])
#loc = open(locfile,'r')
locpot = ChargeDensity(locfile)
#locpot.read_poscar(loc,vel=False)
#locpot.read_chgcar(loc)
den_z = locpot.integrate_z_density()
# Read in lower and upper bounds for average
if len(sys.argv) > 2:
    min = float(sys.argv[2])
    max = float(sys.argv[3])
else:
    print np.average(den_z)
    sys.exit()
# Scale bounds by lenght of z-axis, if not already scaled.
if min > 1.0:
    min = min/locpot.cell_vec[2,2]
if max > 1.0:
    max = max/locpot.cell_vec[2,2]
z_pos = np.linspace(0,1,len(den_z))
# Find indices corresponding to lower and upper bounds
i = 0
while True:
    if z_pos[i] >= min:
        break
    i += 1
j = 0
while True:
    if z_pos[j] >= max:
        break
    j += 1
print i,j,len(z_pos)
print z_pos
print np.average(den_z[min:max])
sys.exit()
