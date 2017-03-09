"""
(updated: 2017/03/09)

Class for planing calculations.
- Simplified formulae adapted from Savitsky 1964, the classic paper.

"""
from __future__ import division
import numpy as np
from astropy import units as u

class Planing(object):
    """
    Compute the lift and drag of a planing plate.
    
    Assumptions:
    -------------
    - Fixed, small angle of attack (small means < 10 deg)
    - Small aspect ratio (< 4)
    - Laminar boundary layer
     - due to the size and slow speed
    
    The inputs should be dimensional.
    """
    
    def __init__(self,tau,lba,b,v,body_m,t_dec):
        ## kinematics
        ### skiing speed
        self.v_sf  = v[0]
        self.v_dec = v[1]
        self.v_ski = v[2]

        ###
        self.t_dec = t_dec

        ## dimensions
        self.lba = lba # aspect ratio
        self.b   = b   # beam (width)
        self.c   = lba * b # chord length
        self.tau = tau.to(u.rad)

        ## body mass
        self.body_m = body_m
        
        ## Constants (water)
        self.rho_water = 1000*u.kg/u.m**3  # density
        self.nu_water  = 1e-6 * u.m**2/u.s # kine. viscosity
        self.g         = 9.8 * u.m/u.s**2

        ####
        self.Re  = (self.v_sf*self.c/self.nu_water).to(u.dimensionless_unscaled)
        self.C_v = (self.v_sf/np.sqrt(self.g*b)).to(u.dimensionless_unscaled)
    
    def drag_coeff(self):
        """
        Drag coefficient, laminar.
        """
        #Re = 1.328 / np.sqrt(self.Re.value) ## laminar
        Re = 0.455 / np.log10(self.Re.value)**2.58 - 1700/self.Re.value
        return Re

    def lift_coeff(self):
        """
        Lift coefficient, planing theory
        """
        pre1 = 0.0120*self.lba**0.5
        pre2 = 0.0055*self.lba**2.5 / self.C_v.value**2
        
        C_L = self.tau.value**1.1 * (pre1 + pre2)
        
        return C_L

    def drag_to_lift_ratio(self):
        """
        Drag to lift force ratio, planing theory
        """
        R = self.tau.value + self.drag_coeff()*self.lba / self.lift_coeff()

        return R
    
    def lift_and_drag_forces(self):
        """
        Returns the lift and drag forces.
        Using the mean speed as an average.
        """
        pref = 0.5*self.rho_water * self.v_sf**2 * self.b**2

        C_L = self.lift_coeff()
        R   = self.drag_to_lift_ratio()

        ## two feet
        F_L = 2 * (pref * C_L).to(u.N)
        F_D = F_L * R
      
        return F_L, F_D

    def drag_power(self):
        """
        Power consumption for overcoming drag
        """
        P_D = (self.lift_and_drag_forces()[1] * self.v_dec).to(u.W)
       
        return P_D

    def decele(self):
        """
        Theoretical deleceration curve.
        From initial velocity (steady flight)
        """
        ## F = m a
        a__ = (self.lift_and_drag_forces()[1] / self.body_m).to(u.cm/u.s**2)
        v_f = self.v_sf - self.t_dec * a__
        
        return v_f

