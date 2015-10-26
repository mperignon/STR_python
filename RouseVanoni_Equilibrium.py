from math import sqrt, log10, pow
import numpy as np

class RouseVanoni_Equilibrium(object):
    """
    RouseVanoni_Equilibrium

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Oct 2015

    ----------------------------------------------

    Computes the Rouse-Vanoni profile of suspended sediment concentration at equilibrium.

    ----------------------------------------------
    Parameters:
    H: Water depth
    eta: Height of the cross-section above the bed
            eta=0 on the bed surface
            eta=H on the water surface
    eta_b: Position near the bed surface where the volumetric concentration of
            suspended sediment is equal to Cb
    u_star: Shear velocity
    vs: Particle settling velocity
    C: Concentration of suspended sediment in the water column at elevation eta
       averaged over turbulence
    Cb: Concentration of suspended sediment in the water column at eta = eta_b
        averaged over turbulence. Cb can be computed with entrainment relationships

    Input:
    b: Non-dimensional distance from the bed surface where Cb is computed
    vs: Settling velocity of the particles (m/s)
    u_star: Shear velocity (m/s)
    spacing: Spacing between points on the vertical profile

    Output:
    Non-dimensional equilibrium profile of suspended sediment concentration.
    The height in the water column is expressed in terms of z = eta/H and the
    volumetric concentration averaged over the turbulence as C/Cb.
        At eta = b, C/Cb = 1
        At eta = 1, C/Cb = 0


    Equations:
    
    Non-dimensional expression of the Rouse-Vanoni profile:
    \frac{C}{Cb} = \left[\frac{(1-z)/z}{(1-b)/b}\right]^{vs/(\kappa*u_\star)}
    where \kappa is VonKarman's constant (0.4), z is a non-dimensional height above the bed
    defined as eta/H, and b is the non-dimensional near-bed distance equal to eta_b/H
    
    
    """
    
    
    def __init__(self, b = 0.05, u_star = 0.2, verbose = True):
        
        self._verbose = verbose
        if verbose:
            print 'Initializing RouseVanoni_Equilibrium...'
            
        self._kappa = 0.4
        
        self._u = float(u_star)
        self._b = float(b)
        self._vs = None
        
        self._z = np.hstack((np.linspace(b,0.95,19), np.array([0.98,0.995,1.0])))
        self._C = np.zeros_like(self._z)
        self._profile = None
        
        
    @property
    def settling_velocity(self):
        """Settling velocity in m/s"""
        return self._vs


    @settling_velocity.setter
    def settling_velocity(self, new_vs):
        """Set the settling velocity.

        Parameters
        ----------
        new_vs : float
            New settling velocity in m/s.
        """
        self._vs = new_vs

    @property
    def RV_profile(self):
        """Rouse-Vanoni Equilibrium profile"""
        
        return self._profile
            
            
    def run(self):
    
        assert self._vs, "Settling velocity must be set before running. Use FallVelocity.py?"
        
        power = ((1-self._z)/self._z) / ((1-self._b)/self._b)
        self._C = np.power(power, self._vs/(self._kappa * self._u))
        
        self._profile = np.vstack([self._z, self._C]).T
            
            
            
    def finalize(self):
                
        print 'The Rouse-Vanoni Equilibrium profile is:'
        print np.flipud(self._profile)
            
            