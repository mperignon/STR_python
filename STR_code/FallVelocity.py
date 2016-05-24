from math import sqrt, log10, pow
import json, os
        

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

    
    def __init__(self, D = 0.1, nu=1e-6, g=9.81, rho_w=1000., rho_s=2650., verbose=False):
        
        self.verbose = verbose
        if self.verbose:
            print 'Initializing FallVelocity...\n'
        
        self.nu = float(nu)
        self.g = float(g)
        self.rho_w = float(rho_w)
        self.rho_s = float(rho_s)
        
        self._grain_size = D / 1000.
        
        if self.verbose:
            print 'Default values set!'
            print '-' * 20
            print 'Grain size:', self.grain_size, 'mm'
            print 'Kinematic viscosity:', self.nu, 'm^2/s'
            print 'Graviatational acceleration:', self.g, 'm/s^2'
            print 'Density of fluid:', self.rho_w, 'Kg/m^3'
            print 'Density of sediment:', self.rho_s, 'Kg/m^3'
        
        # Outputs
        self.vs = None
        self.Re = None
        self.Rf = None
        

    @property
    def grain_size(self):
        """Grain size in mm."""
        return self._grain_size * 1000 # to report in mm


    @grain_size.setter
    def grain_size(self, new_D):
        """Set the grain size.

        Parameters
        ----------
        new_D : float
            New grain size in mm.
        """
        self._grain_size = new_D / 1000.     # convert from mm to m
        
        if self.verbose:
            print 'Grain size set'
            
    
    @property
    def settling_velocity(self):
        """Settling velocity in m/s"""
        
        assert self.vs, "A settling velocity has not been calculated "\
                        "because FallVelocity has not been run."
        return self.vs
    
    @property
    def Reynolds_number(self):
        
        assert self.Re, "A Reynolds number has not been calculated "\
                        "because FallVelocity has not been run."
        
        return self.Re
    
    @property
    def dimensionless_fall_velocity(self):
        
        assert self.Rf, "A dimensionless fall velocity has not been calculated "\
                        "because FallVelocity has not been run."
        
        return self.Rf
        
        
            
    def run(self):
    
        self.R = (self.rho_s - self.rho_w) / self.rho_w
        
        self.calculateReynoldsnumber()     
        self.calculate_dimensionless_fall_velocity()                                       
        self.calculate_settling_velocity()
        
        if self.verbose:
            print 'Running FallVelocity...\n'
            print 'Reynolds Number: %.4g' % self.Re
            print 'Dimensionless fall velocity: %.4g' % self.Rf
            print 'Settling velocity: %.4g m/s' % self.vs
        
        
            
            
    def finalize(self, dir = 'output'):
    
    	if not os.path.exists(dir):
    		os.makedirs(dir)
    
        if self.verbose:
            print 'Finalizing FallVelocity...\n'
            print 'Saved input and output files.'

        input_dict = {
            'Grain_size_mm' : self.grain_size,
            'Kinematic_viscosity' : self.nu,
            'Gravitational_acceleration' : self.g,
            'Density_of_fluid' : self.rho_w,
            'Density_of_sediment' : self.rho_s
        }
        
        with open(os.path.join(dir,'FallVelocity_input.json'), 'w') as f:
            json.dump(input_dict, f, indent=4)
            
        output_dict = {
            'Reynolds_number' : self.Re,
            'Dimensionless_fall_velocity' : self.Rf,
            'Settling_velocity' : self.vs
        }
        
        with open(os.path.join(dir,'FallVelocity_output.json'), 'w') as f:
            json.dump(output_dict, f, indent=4)

                
            
            
            
            

            
    def calculateReynoldsnumber(self):

        self.Re = sqrt(self.R * self.g * self.grain_size) * \
                    self.grain_size / self.nu

        failMessage  = "This equation is only valid for Reynolds numbers "\
                        "below 2.6e+06. The calculated Reynolds number for "\
                        "this particle is %.2g" % self.Re
        assert self.Re <= 2.6e6, failMessage       
  
  
        
    def calculate_dimensionless_fall_velocity(self):
        
        x = log10(self.Re**2)
        y = (-3.76715) + (1.92944*x) - (0.09815*x*x) - \
            (0.00575*x*x*x) + (0.00056*x*x*x*x)
        
        self.Rf = pow(pow(10,y) / self.Re, 1./3)


        
    def calculate_settling_velocity(self):
        
        self.vs = self.Rf * sqrt(self.R * self.g * self.grain_size)
