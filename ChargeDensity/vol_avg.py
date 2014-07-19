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




avg1 = a.avg_density_vol()
avg2 = np.average(a.density)

print "avg1",avg1
print "avg2",avg2



sys.exit()

a.unitcell.scale = 1.0
den_z = a.integrate_z_density()
z_pos = np.linspace(0,a.unitcell.cell_vec[2,2],len(den_z))
macro_z = a.macro_avg_z(p1)
for i in range(len(den_z)):
    print z_pos[i],den_z[i],macro_z[i]

# Calculate bulk and vacuum average, assuming that the bulk is located in the
# 1st half of the cell (along z) and the vacuum is in the second half of the
# cell.
bulk_start = 0.2
bulk_stop = 0.3
vac_start = 0.7
vac_stop = 0.8
bi = int(np.floor(bulk_start*len(den_z)))
bf = int(np.floor(bulk_stop*len(den_z)))
vi = int(np.floor(vac_start*len(den_z)))
vf = int(np.floor(vac_stop*len(den_z)))

bulk_avg = np.average(macro_z[bi:bf])
bulk_std = np.std(macro_z[bi:bf])
#bulk_center = macro_z[int(np.floor(0.25*len(den_z)))]
vac_avg = np.average(macro_z[vi:vf])
vac_std = np.std(macro_z[vi:vf])
#vac_center = macro_z[int(np.floor(0.75*len(den_z)))]

print
print "Bulk_avg_(eV) Bulk_std_(eV) Vac_avg_(eV) Vac_std_(eV)"
print bulk_avg,bulk_std,vac_avg,vac_std
#print "Bulk_avg_(eV) Bulk_center_(eV) Vac_avg_(eV) Vac_center_(eV)"
#print bulk_avg,bulk_center,vac_avg,vac_center
