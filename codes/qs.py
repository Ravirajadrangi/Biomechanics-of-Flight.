"""
(updated: 2017/03/08)

Class for QS model calculations

"""
from __future__ import division
import numpy as np
import pandas as pd
from astropy import units as u

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
        self.wing_ome = (2*np.pi*wing_freq).to(u.s**-1)
        self.wing_amp = wing_amp
        
        ## body kinematics
        self.body_v = body_v
        
        ## precision
        self.Nstep     = Nstep
        self.time_step = np.linspace(0.,1/wing_freq.value,Nstep)*u.s # time steps in one full stroke
        self.Nsegments = Nsegments # how many strips make a wing
        
        ## wing shape : ellipse
        ### 'wing_a' = twice semi-minor axis; 'wing_r' = semi-major axis
        self.wing_a = wing_a
        
        self.wing_r = np.linspace(0.,wing_r.value,Nsegments).reshape(100,1) * u.cm
        self.wing_c = self._ellipse(self.wing_a,wing_r,self.wing_r)
        
        ## Constants
        self.rho_air = 1*u.kg/u.m**3
        
    def angle_of_attack(self):
        """
        Calculate the effective angle of attack variation.
        """
        ##
        wt = self.wing_ome*self.time_step
        
        ## geometric AoA
        ### a "horizontal" wing experiences 90 AoA in downstroke
        alpha_geo = (90 - np.sin(wt.value)*45)*u.deg 
        
        ## correction from body speed
        ### speeds of wing elements
        self.wing_v = (self.wing_amp*(self.wing_r/self.wing_a) * \
                       self.wing_ome*np.sin(wt.value)).to(u.cm/u.s)
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
        ### note that both terms multiple a sine (corr ang defined w.r.t. v_wing)
        F_vert = pref*self.v_tot**2*C_L*np.sin(np.deg2rad(self.alpha_cor))
        F_hori = pref*self.v_tot**2*C_D*np.sin(np.deg2rad(self.alpha_cor))
        
        ## integration; two wings
        F_V = 2 * np.trapz(F_vert,self.wing_r,axis=0).to(u.N)
        F_H = 2 * np.trapz(F_hori,self.wing_r,axis=0).to(u.N)

        ## upstroke = 0
        ### but storke upstroke
        F_V_,F_H_ = np.copy(F_V), np.copy(F_H)
        
        F_V[self.Nstep/2:] = 0
        F_H[self.Nstep/2:] = 0
    
        return F_V, F_H, F_V_, F_H_

    def transl_d_power(self):
        """
        Power consumption for overcoming drag
        """
        ## same as force
        pref = 0.5*self.rho_air*self.wing_c
        C_L,C_D = self.lift_drag_coeff()

        self.v_tot = np.sqrt(self.wing_v**2 + self.body_v**2)

        ## power = F_D * v_tot
        P_drag = pref*self.v_tot**2*C_D * self.v_tot

        ## integration; two wings
        P_D = 2 * np.trapz(P_drag,self.wing_r,axis=0).to(u.W)

        ## upstroke = 0
        ### but keeps upstroke
        P_D_ = np.copy(P_D)
        
        P_D[self.Nstep/2:] = 0

        return P_D, P_D_
        
        
    #########################
    # Supplementary methods #
    def _ellipse(self,a2,b,r):
        """
        Given radius returns chord length
        """
        return a2 * np.sqrt(1 - (r/b)**2)

