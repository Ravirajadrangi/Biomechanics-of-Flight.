"""
(updated: 2017/03/09)

Class for brake calculations. Mostly modified from qs.py

"""
from __future__ import division
import numpy as np
from astropy import units as u

class Brake(object):
    """
    Compute the horizontal drag force and power in steady braking.
    
    Assumptions:
    -------------
    - Fixed angle of attack. Use the body orientation angle.
    - Wings are half-ellipses
     - 2a = chord length at shoulder
     -  b = radius (span) of the wing
    - Non-flapping wings
     - so u_e = u_body
    - Constant deceleration
     - i.e. u_body(t) = u_body(0) - a t
    
    The inputs should be dimensional.
    """
    
    def __init__(self,wing_a,wing_r,alpha,body_mv,body_sv,t_dec,body_m):
        ## kinematics
        ### 'mean' body speed
        self.body_mv = body_mv

        ### 'initial' body speed = steady flight speed
        self.body_sv = body_sv

        ### time to decelerate to mean speed (obs.)
        self.t_dec = t_dec
               
        ## wing
        ### 'wing_a' = twice semi-minor axis; 'wing_r' = semi-major axis
        self.wing_a = wing_a
        self.wing_r = wing_r

        ### wing area = \pi a b / 2 = \pi 2a b / 4
        self.wing_area = np.pi*wing_a*wing_r/4

        ### AoA
        self.alpha = alpha

        ## body mass
        self.body_m = body_m
        
        ## Constants
        self.rho_air = 1*u.kg/u.m**3
    
    def drag_coeff(self):
        """
        Drag coefficient, at the given AoA.
        """
        ## alpha always > 0
        a_D = 0.0073*self.alpha.value+3.1416
        C_D =  8.3171 + 8.1909 * np.cos(a_D)
                
        return C_D
    
    def transl_force(self):
        """
        Returns the translational drag force for the mean speed
        Using the initial speed
        """
        pref = 0.5*self.rho_air*self.wing_area
        C_D = self.drag_coeff()
        F_D = 2 * (pref*self.body_sv**2*C_D).to(u.N)
      
        return F_D

    def transl_d_power(self):
        """
        Power consumption for overcoming drag
        """
        P_D = (self.transl_force() * self.body_mv).to(u.W)
       
        return P_D

    def decele(self):
        """
        Theoretical deleceration curve.
        From initial velocity (steady flight)
        """
        ## F = m a
        a__ = (self.transl_force() / self.body_m).to(u.cm/u.s**2)
        v_f = self.body_sv - self.t_dec * a__
        
        return v_f
