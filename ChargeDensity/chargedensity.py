#!/usr/bin/env python

# ChargeDensity Class v1.0 07-19-2014 Jeff Doak jeff.w.doak@gmail.com

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
        """
        ChargeDensity(input_,format_)
            input_ - Can be a ChargeDensity object, a string representing a file
            name containing a CHGCAR or LOCPOT formatted file, a file-type
            object containing density data in the format of a CHGCAR file, or
            None.
            format_ - can be None, locpot, or chgcar. Determines whether or not
            to divide the data in the file by the volume of the unit cell. None
            and locpot do not divide the data, while chgcar does.
        """
        self.unitcell = UnitCell()
        if isinstance(input_,ChargeDensity):
            # Return a new ChargeDensity object that is a copy of input_
            self.unitcell = UnitCell(input_.unitcell)
            self.density = input_.density
        elif isinstance(input_,str):
            try:
                input_ = open(input_,'r')
            except IOError:
                print "Error reading input file."
                print "Empty ChargeDensity will be returned."
                input_ = None
        #if isinstance(input_,file):
        if input_ is not None:
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
        if input_ is None:
            # Create empty chargedensity.
            self.density = np.zeros((1,1,1))

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

    def periodize_density(self):
        """
        Returns a copy of the density with increased size (nx+1,ny+1,nz+1) with
        density[:,:,-1] = density[:,:,0], density[:,-1,:] = density[:,0,:],
        and density[-1,:,:] density[0,:,:], making the density function periodic
        in and of itself. This is important for integrating the density over the
        volume of the unit cell.
        """
        (nx,ny,nz) = np.shape(self.density)
        per_density = np.zeros((nx+1,ny+1,nz+1))
        per_density[0:nx,0:ny,0:nz] = self.density
        per_density[-1,:,:] = per_density[0,:,:]
        per_density[:,-1,:] = per_density[:,0,:]
        per_density[:,:,-1] = per_density[:,:,0]
        return per_density

    def periodize_den_xy(self):
        """
        Returns a copy of the density with increased size (nx+1,ny+1,nz) with
        density[-1,:,:] = density[0,:,:], and density[:,-1,:] =  density[:,0,:],
        making the density function periodic in and of itself along x and y.
        This is important for integrating the density in the x-y plane.
        """
        (nx,ny,nz) = np.shape(self.density)
        per_density = np.zeros((nx+1,ny+1,nz))
        per_density[0:nx,0:ny,0:nz] = self.density
        per_density[-1,:,:] = per_density[0,:,:]
        per_density[:,-1,:] = per_density[:,0,:]
        return per_density

    def integrate_z_density(self):
        """
        Returns a numpy array containing the density integrated in the xy-plane
        as a function of position along the z-axis. Uses the trapezoid rule to
        integrate the density.
        """
        from scipy.integrate import simps
        den = self.periodize_den_xy()
        (nx,ny,nz) = np.shape(den)
        a = np.linalg.norm(self.unitcell.cell_vec[0])
        b = np.linalg.norm(self.unitcell.cell_vec[1])
        den_yz = simps(den,dx=float(a/(nx-1)),axis=0)/a
        den_z = simps(den_yz,dx=float(b/(ny-1)),axis=0)/b
        return den_z

    def avg_density_vol(self):
        """
        Averages the density over the entire volume of the calculation cell.
        """
        from scipy.integrate import trapz
        den = self.periodize_density()
        (nx,ny,nz) = np.shape(den)
        self.unitcell.scale = 1.0
        a = np.linalg.norm(self.unitcell.cell_vec[0])
        b = np.linalg.norm(self.unitcell.cell_vec[1])
        c = np.linalg.norm(self.unitcell.cell_vec[2])
        vol = np.linalg.det(self.unitcell.cell_vec)
        jac = self.unitcell.cell_vec*np.outer(np.array((1./a,1./b,1./c)),np.ones(3))
        jac_det = np.linalg.det(jac)
        den_yz = trapz(den,dx=float(1./(nx-1)),axis=0)
        den_z = trapz(den_yz,dx=float(1./(ny-1)),axis=0)
        avg_den = trapz(den_z,dx=float(1./(nz-1)),axis=0)
        avg_den = avg_den*jac_det/vol
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

    def interpolate(self,point):
        """
        Finds the density at an arbitrary point in direct coordinates by
        tri-linear interpolation of the surrounding density values.
        """
        # From Scott Kirklin
        # Make sure the point is located within the interior of the unit cell.
        point = [ p%1 for p in point ]
        point = [ p%1 for p in point ]
        # Find lower bounding points x0,y0,z0
        x0,y0,z0 = (int(np.floor(point[0]*self.mesh[0])),
                    int(np.floor(point[1]*self.mesh[1])),
                    int(np.floor(point[2]*self.mesh[2])))
        # Find upper bounding points x1,y1,z1
        x1,y1,z1 = ((x0+1)%self.mesh[0],
                    (y0+1)%self.mesh[1],
                    (z0+1)%self.mesh[2])
        # We now have the bouding points: (x0,y0,z0), (x1,y0,z0), (x0,y1,z0),
        # (x0,y0,z1), (x1,y1,z0), (x1,y0,z1), (x0,y1,z1), (x1,y1,z1) that form a
        # box surrounding point. Now find the fractional position of point in
        # the boundiing box.
        x,y,z   =  ((point[0]*self.mesh[0])%1,
                    (point[1]*self.mesh[1])%1,
                    (point[2]*self.mesh[2])%1)
        # The interpolated value is a linear combination of the values at each
        # corner of the box.
        interp_val = (self.density[x0,y0,z0]*(1-x)*(1-y)*(1-z) +
                      self.density[x1,y0,z0]*x*(1-y)*(1-z) +
                      self.density[x0,y1,z0]*(1-x)*y*(1-z) +
                      self.density[x0,y0,z1]*(1-x)*(1-y)*z +
                      self.density[x1,y1,z0]*x*y*(1-z) +
                      self.density[x1,y0,z1]*x*(1-y)*z +
                      self.density[x0,y1,z1]*(1-x)*y*z +
                      self.density[x1,y1,z1]*x*y*z)
        return interp_val

    def dens_diff_z(self,dens1,dens2):
        """
        Takes the difference between two 1-d densities using the grid spacing
        corresponding to dens1. Assumes dens1 and dens2 are splines.
        """
        z_pos = np.linspace(0,self.unitcell.cell_vec[2,2],200)
        diff = dens2(z_pos) - dens1(z_pos)
        return z_pos,diff

    def spherical_average(self,center,r_max,n_pts=10):
        """
        Averages the density in a sphere around the point center. Returns an
        array containing the spherically averaged density over a series of
        shells with increasing radii centered at the point center.
        """
        # Physics notation for spherical coordinates are use:
        # rho - radial distance
        # theta - inclination
        # phi - azimuthal angle
        from scipy.integrate import dblquad
        def func(theta,phi,rho,center,dens):
            """
            Function to find the value of the density at a point equal to
            center+Vector(rho).
            """
            # Determine direct rectilinear coordinates relative to center from
            # spherical coordinates away from center.
            rho_x = rho*np.cos(phi)*np.sin(theta)
            rho_y = rho*np.sin(phi)*np.sin(theta)
            rho_z = rho*np.cos(theta)
            # Determine absolute direct coordinates from relative coordinates
            # and center.
            r_x = rho_x + center[0]
            r_y = rho_y + center[1]
            r_z = rho_z + center[2]
            # Find the value of the density at the point r.
            val = dens.interpolate([r_x,r_y,r_z])
            # Weight the value by the spherical coordinate Jacobian.
            val = val*rho*rho*np.sin(theta)
            return val

        r_list = np.linspace(0,r_max,n_pts)
        val_list = np.zeros_like(r_list)
        for i in range(len(r_list)):
            rho = r_list[i]
            integral,error = dblquad(func,0,2*np.pi,lambda g:0, lambda h: np.pi,
                                     args=(rho,center,self))
            val_list[i] = integral
        return r_list,val_list

    def spherical_vol_avg(self,center,r_max):
        """
        Averages the density over the volume of a sphere centered at the point
        center, with a radis r_max.
        """
        # Physics notation for spherical coordinates are use:
        # rho - radial distance
        # theta - inclination
        # phi - azimuthal angle
        from scipy.integrate import tplquad
        def func(theta,phi,rho,center,dens):
            """
            Function to find the value of the density at a point equal to
            center+Vector(rho).
            """
            # Determine direct rectilinear coordinates relative to center from
            # spherical coordinates away from center.
            rho_x = rho*np.cos(phi)*np.sin(theta)
            rho_y = rho*np.sin(phi)*np.sin(theta)
            rho_z = rho*np.cos(theta)
            # Determine absolute direct coordinates from relative coordinates
            # and center.
            r_x = rho_x + center[0]
            r_y = rho_y + center[1]
            r_z = rho_z + center[2]
            # Find the value of the density at the point r.
            val = dens.interpolate([r_x,r_y,r_z])
            # Weight the value by the spherical coordinate Jacobian.
            val = val*rho*rho*np.sin(theta)
            return val

        integral,error = tplquad(func,0,r_max,lambda a:0,lambda b:2*np.pi,
                                 lambda c,d:0,lambda e,f: np.pi,args=(center,self))
        return integral


if __name__ == "__main__":
    a = ChargeDensity("LOCPOT")

