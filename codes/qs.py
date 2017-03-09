"""
(updated: 2017/03/08)

Class for QS model calculations

"""

from __future__ import division
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

class QS(object):
    """
    Compute the horizontal, vertical, and total forces from QS model (two wings).
    Drag power is calculated too.
    
    Assumptions:
    -------------
    - Sinusoidal wing kinematics (both strokes and rotations)
     - The maximum geometric angle of attack is 45 deg
    - Wings are half-ellipses
     - 2a = chord length at shoulder
     -  b = radius (span) of the wing
    - Constant motions (including body speed)
    
    The inputs should be dimensional.
    """
    
    def __init__(self,wing_freq,wing_amp,wing_a,wing_r,body_v,
                 Nstep=200,Nsegments=100):
        ## wing kinematics
        self.wing_ome = 2*np.pi*wing_freq
        self.wing_amp = wing_amp
        
        ## body kinematics
        self.body_v = body_v
        
        ## precision
        self.time_step = np.linspace(0.,1/wing_freq.value,Nstep) # time steps in one full stroke
        self.Nsegments = Nsegments # how many strips make a wing
        
        ## wing shape : ellipse
        ### 'wing_a' = twice semi-minor axis; 'wing_r' = semi-major axis
        self.wing_r = np.linspace(0.,wing_r.value,Nsegments).reshape(100,1) * u.cm
        self.wing_c = self._ellipse(wing_a,wing_r,self.wing_r)
        
        ## Constants
        self.rho_air = 1*u.kg/u.m**3
        
    def angle_of_attack(self):
        """
        Calculate the effective angle of attack variation.
        """
        ## geometric AoA
        ### a "horizontal" wing experiences 90 AoA in downstroke
        alpha_geo = (90 - np.sin(self.wing_ome.value*self.time_step)*45)*u.deg 
        
        ## correction from body speed
        ### speeds of wing elements
        self.wing_v = (self.wing_ome * self.wing_r).to(u.cm/u.s)
        self.alpha_cor = np.mod(90-np.rad2deg(np.arctan(self.wing_v.value \
                                                       /self.body_v.value)),180)*u.deg
        
        return alpha_geo - self.alpha_cor
    
    def lift_drag_coeff(self):
        """
        Lift and drag coefficients, at a given AoA.
        
        Note: C_L are scaled up by 2x from hummingbird to pigeon.
        """
        ## initialize
        alpha = self.angle_of_attack()
        C_L,C_D = np.zeros(alpha.shape),np.zeros(alpha.shape)
        
        for ((i,j),alp) in np.ndenumerate(alpha):
            if alp >= 0:
                a_L = 0.0301*alp+4.7124
                a_D = 0.0073*alp+3.1416
                C_L[i,j] = (0.0031 + 1.5842 * np.cos(a_L)) * 2
                C_D[i,j] =  8.3171 + 8.1909 * np.cos(a_D)
            else:
                a_L = 0.0332*alp+4.6963
                a_D = 0.0281*alp+3.1277
                C_L[i,j] = (0.0028 + 1.1251 * np.cos(a_L)) * 2
                C_D[i,j] =  1.1993 + 1.0938 * np.cos(a_D)
        
        return C_L, C_D
    
    def transl_forces(self):
        """
        Returns the translational forces at each time stamp
        """
        pref = 0.5*self.rho_air*self.wing_c
        C_L,C_D = self.lift_drag_coeff()
        
        ## total wing speed relative to air
        self.v_tot = np.sqrt(self.wing_v**2 + self.body_v**2) 
        
        ## transform each element from lift/drag to vert/hori
        ### note that both terms multiple a cosine of 
        F_vert = pref*self.v_tot**2*C_L*np.abs(np.cos(np.deg2rad(self.alpha_cor)))
        F_hori = pref*self.v_tot**2*C_D*np.abs(np.cos(np.deg2rad(self.alpha_cor)))
        
        ## integration
        F_L = np.trapz(F_vert,self.wing_r,axis=0).to(u.N)
        F_D = np.trapz(F_hori,self.wing_r,axis=0).to(u.N)
    
        return F_L, F_D
        
        
    #########################
    # Supplementary methods #
    def _ellipse(self,a2,b,r):
        """
        Given radius returns chord length
        """
        return a2 * np.sqrt(1 - (r/b)**2)
