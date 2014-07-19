#!/usr/bin/env python

# ChargeDensity Class v0.1 09-19-2012 Jeff Doak jeff.w.doak@gmail.com

# Changelog:

from unitcell import *
import numpy as np
import sys

class ChargeDensity():
    """
    Class to read in, store, and manipulate charge density files.
    
    Instance attributes:
    self.unitcell - unitcell object of the POSCAR corresponding to the density
    self.density  - density (charge or electrostatic potential) as a function of
                    position in the unit cell
    """

    def __init__(self,input_=None,format_=None):
        self.unitcell = UnitCell()
        if isinstance(input_,str):
            try:
                input_ = open(input_,'r')
            except IOError:
                print "Error reading input file."
                print "Empty ChargeDensity will be returned."
                input_ = None
        if isinstance(input_,file):
            if format_ == None or format_ == "locpot":
                self.unitcell.read_poscar(input_,vel=False)
                self.read_chgcar(input_)
            elif format_ == "chgcar":
                self.unitcell.read_poscar(input_,vel=False)
                self.read_chgcar(input_)
                self.scale_density()
            else:
                print "Unknown format."
                print "Empty ChargeDensity will be returned."
                input_ = None
        elif isinstance(input_,ChargeDensity):
            # Return a new ChargeDensity object that is a copy of input_
            self.unitcell = UnitCell(input_.unitcell)
            self.density = input_.density

    def read_chgcar(self,input_):
        """
        Reads in a charge-density from file input_. Assumes that the atomic
        positions have been read from the file already.
        """
        line = input_.readline().split()
        ngx = int(line[0]); ngy = int(line[1]); ngz = int(line[2])
        density = np.zeros((ngx,ngy,ngz))
        x = 0; y = 0; z = 0
        for line in input_:
            line = line.split()
            for i in range(len(line)):
                density[x,y,z] = float(line[i])
                x += 1
                if x >= ngx:
                    x = 0
                    y += 1
                    if y >= ngy:
                        y = 0
                        z += 1
                        if z >= ngz:
                            self.density = density
                            return
        self.density = density
        return

    def scale_density(self,scale=None):
        """
        Scales the electron density by a constant. If scale is None, the inverse
        of the volume of the unit cell is used as the scale.
        """
        if scale == None:
            self.scale = 1.0
            scale = 1.0/np.linalg.det(self.unitcell.cell_vec)
        self.density = self.density*scale

    def integrate_z_density(self):
        """
        Returns a numpy array containing the density integrated in the xy-plane
        as a function of position along the z-axis. Uses the trapezoid rule to
        integrate the density.
        """
        from scipy.integrate import trapz
        ngx = len(self.density)
        ngy = len(self.density[0])
        ngz = len(self.density[0,0])
        a = np.linalg.norm(self.unitcell.cell_vec[0])
        b = np.linalg.norm(self.unitcell.cell_vec[1])
        den_yz = trapz(self.density,dx=float(a/ngx),axis=0)/a
        den_z = trapz(den_yz,dx=float(b/ngy),axis=0)/b
        return den_z

    def avg_density_vol(self):
        """
        Averages the density over the entire volume of the calculation cell.
        """
        from scipy.integrate import trapz
        ngx = len(self.density)
        ngy = len(self.density[0])
        ngz = len(self.density[0,0])
        a = np.linalg.norm(self.unitcell.cell_vec[0])
        b = np.linalg.norm(self.unitcell.cell_vec[1])
        c = np.linalg.norm(self.unitcell.cell_vec[2])
        vol = np.linalg.det(self.unitcell.cell_vec)
        jac = self.unitcell.cell_vec*np.outer(np.array((1./a,1./b,1./c)),np.ones(3))
        jac_det = np.linalg.det(jac)
        den_yz = trapz(self.density,dx=float(1./ngx),axis=0)
        den_z = trapz(den_yz,dx=float(1./ngy),axis=0)
        avg_den = trapz(den_z,dx=float(1./ngz),axis=0)
        avg_den = avg_den/jac_det/vol
        return avg_den

    def macro_avg_z(self,lat_const=None):
        """
        Returns a numpy array containing the macroscopic average of a density
        along the z-axis that has been integrated in the xy-plane. If no lattice
        constant to average over is provided, the length of the smallest unit
        cell vector will be used.
        """
        from numpy.fft import rfft,irfft
        if lat_const == None:
            lat_const = 1.0
        # Convert periods to direct units, if given in cartesian
        if lat_const > 1.0:
            lat_const = lat_const/self.unitcell.cell_vec[2,2]
        den_z = self.integrate_z_density()
        z_pos = np.linspace(0,1,len(den_z))
        # Find index of lowest lower bound
        low = 1.-lat_const/2.
        i = len(z_pos)-1
        while True:
            if z_pos[i] <= low:
                i += 1
                break
            i -= 1
        #Find index of lowest upper bound
        high = lat_const/2.
        j = 0
        while True:
            if z_pos[j] >= high:
                j -= 1
                break
            j += 1
        # Compare with convolution via fourier transforms
        s = np.zeros(len(den_z))
        s[0:j+1] = np.ones(j+1)
        s[i:] = np.ones(len(s[i:]))
        s = s/float(np.sum(s))
        f_z = rfft(den_z)
        f_s = rfft(s)
        macro_z = irfft(f_z*f_s)
        return macro_z

    def macro_mismatch(self,p1,p2):
        """
        Performs double convolution with two different periods to calculate
        macroscopic average of a charge density along the z-axis.
        """
        from numpy.fft import rfft,irfft
        # Convert periods to direct units, if given in cartesian
        if p1 > 1.0:
            p1 = p1/self.unitcell.cell_vec[2,2]
        if p2 > 1.0:
            p2 = p2/self.unitcell.cell_vec[2,2]
        # Create xy-plane averaged density
        micro_z = self.integrate_z_density()
        # Create convolutions
        z_pos = np.linspace(0,1,len(micro_z))
        # Find index of lowest lower bound for p1
        low = 1.-p1/2.
        i1 = len(z_pos)-1
        while True:
            if z_pos[i1] <= low:
                i1 += 1
                break
            i1 -= 1
        #Find index of lowest upper bound for p1
        high = p1/2.
        j1 = 0
        while True:
            if z_pos[j1] >= high:
                j1 -= 1
                break
            j1 += 1
        # Find index of lowest lower bound for p2
        low = 1.-p2/2.
        i2 = len(z_pos)-1
        while True:
            if z_pos[i2] <= low:
                i2 += 1
                break
            i2 -= 1
        #Find index of lowest upper bound for p2
        high = p2/2.
        j2 = 0
        while True:
            if z_pos[j2] >= high:
                j2 -= 1
                break
            j2 += 1
        conv1 = np.zeros(len(micro_z))
        conv1[0:j1+1] = np.ones(j1+1)
        conv1[i1:] = np.ones(len(conv1[i1:]))
        conv2 = np.zeros(len(micro_z))
        conv2[0:j2+1] = np.ones(j2+1)
        conv2[i2:] = np.ones(len(conv2[i2:]))
        # Perform convolutions in Fourier Space
        f_micro_z = rfft(micro_z)
        f_conv1 = rfft(conv1)
        f_conv2 = rfft(conv2)
        f_macro = f_conv2*f_conv1*f_micro_z
        macro_z = irfft(f_macro)/float(np.sum(conv1))/float(np.sum(conv2))
        return macro_z

    def spline_density_z(self,dens):
        """
        Fits a 1-d density to a cubic spline. Useful for taking
        difference between densities on different grids.
        """
        from scipy.interpolate import UnivariateSpline
        z_pos = np.linspace(0,self.unitcell.cell_vec[2,2],len(dens))
        spline_dens = UnivariateSpline(z_pos,dens,s=0)
        return spline_dens

    def dens_diff_z(self,dens1,dens2):
        """
        Takes the difference between two 1-d densities using the grid spacing
        corresponding to dens1. Assumes dens1 and dens2 are splines.
        """
        z_pos = np.linspace(0,self.unitcell.cell_vec[2,2],200)
        diff = dens2(z_pos) - dens1(z_pos)
        return z_pos,diff

if __name__ == "__main__":
    a = ChargeDensity("LOCPOT")

