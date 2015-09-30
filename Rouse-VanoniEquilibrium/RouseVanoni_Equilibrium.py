from math import sqrt, log10, pow
import numpy as np

class RouseVanoni_Equilibrium:
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
    
    
    def __init__(self, b = None, u_star = None, vs = None, verbose = True):
        
        if verbose:
            print 'Initializing RouseVanoni_Equilibrium...'
            
        self.kappa = 0.4
        self.verbose = verbose
        
        self.u = float(u_star)
        self.b = float(b)
        self.vs = float(vs)
        
        self.z = np.hstack((np.linspace(b,0.95,19), np.array([0.98,0.995,1.0])))
        self.C = np.zeros_like(self.z)

        self.run()
        self.finalize()
            
            
            
            
            
    def run(self):
        
        power = ((1-self.z)/self.z) / ((1-self.b)/self.b)
        self.C = np.power(power, self.vs/(100 * self.kappa * self.u))
            
            
            
            
    def finalize(self):
    
        if self.verbose:
            string = ''
            for i in range(-1,-len(f.z)-1,-1):
                string += '[%.3f, %.3f]\n' % (f.z[i], f.C[i]) 
                
            print 'The Rouse-Vanoni Equilibrium profile is:'
            print string
            
            