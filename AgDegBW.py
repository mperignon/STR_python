from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

import shared_tools as tools

class AgDegBW(object):
    
    def __init__(self,
                 friction = "Manning",
                 Qf = 400,
                 If = 0.1,
                 S = 0.0002,
                 initdepth = 15,
                 B = 60,
                 L = 50000,
                 dx = 833,
                 Gtf = 187059,
                 D90 = 0.6,
                 porosity = 0.4,
                 kc = 75,
                 phi = 1,
                 alpha_r = 8.1,
                 alpha_t = 8.,
                 tau_c = 0.047,
                 n = 1.5,
                 Cz = 22,
                 alpha_u = 1,
                 dt = 0.1,
                 totalT = 120,
                 g = 9.81,
                 rho_w = 1000,
                 rho_s = 2650):
    
    
        self._qw = float(Qf) / float(B)

        self._friction = friction
        self._initdepth = initdepth

        self._Qf = float(Qf) # characteristic flood discharge
        self._If = float(If) # flood intermittency
        self._B = float(B) # channel width

        self._porosity = float(porosity)
        self._S = float(S)
        self._Gtf =  float(Gtf) # imposed annual sediment transport into reach
        self._L = float(L) # length of fluvial reach
        self._dx = float(dx)
        
        self._D90 = float(D90)/1000
        self._kc = float(kc)/1000

        self._nodes = int(self._L/self._dx) +1 # number of sub-reaches

        self._dt = float(dt) # years
        self._totalT = float(totalT)

        self._alpha_u = float(alpha_u) # upwind coefficient [1=full upwind, 0.5=central difference]
        self._phi = float(phi) # boundary shear stress due to skin friction (phi <= 1)
        self._alpha_r = float(alpha_r) # coefficient in Manning resistance relation

        # parameters for load relation (Wong and Parker, 2006)
        self._alpha_t = float(alpha_t) # coefficient in sed transport relation
        self._tau_c = float(tau_c)
        self._n = float(n) # exponent in sed transport relation

        # Chezy
        self._Cf = 1/float(Cz)**2

        # Constants
        self._g = float(g)
        self._sec_to_year = 60 * 60 * 24 * 365.25

        # Derivatives
        self._R = (float(rho_s) - float(rho_w)) / float(rho_w)

        # arrays
        self._x, self._eta, self._H, self._ksi = tools.create_surface_arrays(self._nodes,
                                                                                 self._dx,
                                                                                 self._S)

        self._qb = np.zeros_like(self._x)
        self._tau = np.zeros_like(self._x)
        self._Ht = np.zeros_like(self._x)
        
        
    @property
    def bed_elevation(self):
        return self._eta
    
    @property
    def water_surface_elevation(self):
        return self._ksi
    
    @water_surface_elevation.setter
    def water_surface_elevation(self, new_ksi):
        self._ksi = new_ski
        
    @property
    def dam_height(self):
        return self._initdepth
    
    @property
    def slope(self):
        return self._Sl
    
    @slope.setter
    def slope(self, new_slope):
        self._Sl = new_slope
        
    @property
    def bed_shear_stress(self):
        return self._tau
    
    @bed_shear_stress.setter
    def bed_shear_stress(self, new_tau):
        self._tau = new_tau
        
    @property
    def distance(self):
        return self._x
        

    def Backwater_calculator(self):
        
        self._H[-1] = self._initdepth - self._eta[-1]
        
        for i in range(-1,-self._nodes,-1):
        
            Fr2p = self._qw**2 / (self._g * self._H[i]**3)

            if self._friction == "Manning":
                Cf = (1/self._alpha_r**2) * (self._kc / self._H[i])**(1./3)

            if self._friction == "Chezy":
                Cf = self._Cf

            fnp = (self._eta[i-1] - self._eta[i] - (Cf * Fr2p * self._dx)) / (1 - Fr2p)
            
            Hp = self._H[i] - fnp
            
            Fr2 = self._qw**2 / (self._g * Hp**3)
            
            fn = (self._eta[i-1] - self._eta[i] - (Cf * Fr2 * self._dx)) / (1 - Fr2)
            
            self._H[i-1] = self._H[i] - 0.5*(fnp + fn)
            
            
        self._ksi = self._eta + self._H
            


    def advance_one_timestep(self):

        
        if self._friction == "Manning":
            Cf = (1/self._alpha_r**2) * (self._kc / self._H)**(1./3)

        if self._friction == "Chezy":
            Cf = self._Cf
        
        self._tau = Cf * (self._qw / self._H)**2 / (self._g * self._R * self._D90)
        

        self._load = self._phi * self._tau - self._tau_c
        self._load[self._load <= 0] = 0
        self._load[self._load > 0] = self._alpha_t * self._load[self._load > 0]**self._n


        self._qb = self._load * self._D90 * sqrt(self._R * self._g * self._D90)

        _qtf = self._Gtf / (self._B * (1+self._R) * self._If * self._sec_to_year)

        self._dq = tools.central_difference(self._qb,
                                            self._dx,
                                            alpha = self._alpha_u,
                                            left_boundary = _qtf)
        
        
        self._eta = self._eta + \
            (self._dt * self._sec_to_year / (1-self._porosity)) * self._If * self._dq




    def run(self):
        
        self._ksi[-1] = self._initdepth
        self.Backwater_calculator()
        
        for t in np.arange(0, self._totalT + self._dt, self._dt):
            
            self.advance_one_timestep()
            self.Backwater_calculator()
            

            if (t%1)==0:
                print t
                
                plt.plot(self._x, self._ksi, self._x, self._eta)
                plt.show()