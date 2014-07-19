#!/usr/bin/python

# atomic_locpot.py v0.1 01/180/2013 Jeff Doak jeff.w.doak@gmail.com

# This script finds the average electrostatic potential over the volume of a
# sphere centered around each atom in the POSCAR file.

from chargedensity import *
import numpy as np
import sys

a = ChargeDensity("LOCPOT")

a.unitcell.convention = "Direct"
a.unitcell.scale = 1.0

radius = float(sys.argv[1]) # radius of volume in angstroms
radius = radius/a.unitcell.cell_vec[0,0]  # This is only valid for cubic cells!

# Set atom names if provided as input argument
if len(sys.argv) > 2:
    names = []
    for i in range(len(a.unitcell.atom_types)):
        names.append(str(sys.argv[2+i]))
    a.unitcell.atom_names = names

#for i in range(a.unitcell.num_atoms):
for i in range(1):
    center = a.unitcell.atom_positions[i]
    V_avg = a.spherical_vol_avg(center,radius)
    name = a.unitcell.atom_names[i]
    x = a.unitcell.atom_positions[i][0]
    y = a.unitcell.atom_positions[i][1]
    z = a.unitcell.atom_positions[i][2]
    print name,x,y,z,V_avg
sys.exit()
