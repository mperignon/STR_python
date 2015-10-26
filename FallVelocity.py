from math import sqrt, log10, pow
        

class FallVelocity(object):
    """
    FallVelocity

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Oct 2015

    ----------------------------------------------

    Computes the settling velocity, vs, of a particle with the formulation of Dietrich (1982).

    Only valid for Reynold's numbers less than or equal to 2.5 x 10**6. If Re is greater than
    this upper limit, warns and exits.

    Dietrich, E. W., 1982, Settling velocity of natural particles,
    Water Resources Research, 18(6), 1626-1982.

    ----------------------------------------------
    Input:
    D: Sediment particle size (mm)
    nu: Kinematic viscosity of the liquid (m**2/s)
    g: Acceleration due to gravity (m/s**2)
    rho_w: Density of water (Kg/m**3)
    rho_s: Density of sediment(Kg/m**3)

    Output:
    vs: Particle settling velocity (m/s)
    Re: Particle Reynolds number
    Rf: Dimensionless fall velocity of the particle


    Equations:
    R = (rho_s - rho_w) / rho_w

    Re = sqrt(R*g*D)*D / nu

    Rf = vs / sqrt(R*g*D)
    """

    
    def __init__(self, nu=1e-6, g=9.81, rho_w=1000., rho_s=2650., verbose=True):
        
        self._verbose = verbose
        if self._verbose:
            print 'Initializing FallVelocity...'
        
        self._nu = nu
        self._g = g
        self._rho_w = float(rho_w)
        self._rho_s = float(rho_s)
        
        self._R = (self._rho_s - self._rho_w) / self._rho_w
        
        self._grain_size = None
        self._vs = None
        

    @property
    def grain_size(self):
        """Grain size in mm."""
        return self._grain_size


    @grain_size.setter
    def grain_size(self, new_D):
        """Set the grain size.

        Parameters
        ----------
        new_D : float
            New grain size in mm.
        """
        self._grain_size = new_D / 1000     # convert from mm to m
            
    
    @property
    def settling_velocity(self):
        """Settling velocity in m/s"""
        
        assert self._vs, "A settling velocity has not been calculated "\
                        "because FallVelocity has not been run."
        return self._vs
        
        
            
    def run(self):
    
        assert self._grain_size, "Grain size must be set before running."
        
        self.calculate_Reynolds_number()     
        self.calculate_dimensionless_fall_velocity()                                       
        self.calculate_settling_velocity()
        
        
            
            
    def finalize(self):
    
        if self._verbose:
            print 'Reynolds Number: %.4g' % self._Re
            print 'Dimensionless fall velocity: %.4g' % self._Rf
            print 'Settling velocity: %.4g m/s' % self._vs

            
    def calculate_Reynolds_number(self):

        self._Re = sqrt(self._R * self._g * self._grain_size) * \
                    self._grain_size / self._nu

        failMessage  = "This equation is only valid for Reynolds numbers "\
                        "below 2.6e+06. The calculated Reynolds number for "\
                        "this particle is %.2g" % self._Re
        assert self._Re <= 2.6e6, failMessage       
  
  
        
    def calculate_dimensionless_fall_velocity(self):
        
        x = log10(self._Re**2)
        y = (-3.76715) + (1.92944*x) - (0.09815*x*x) - \
            (0.00575*x*x*x) + (0.00056*x*x*x*x)
        
        self._Rf = pow(pow(10,y) / self._Re, 0.33)
        
   
        
    def calculate_settling_velocity(self):
        
        self._vs = self._Rf * sqrt(self._R * self._g * self._grain_size)
