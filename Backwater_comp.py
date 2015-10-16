import numpy as np
from math import sqrt

import shared_tools as tools

class Backwater(object):
    
    def __init__(self,
                 friction = "Chezy",
                 qw = 5.7,
                 S = 0.00025,
                 initdepth = 30,
                 L = 123000,
                 dx = 1000,
                 alpha = 8.1,
                 n = 2,
                 D90 = .91,
                 D50 = .27,
                 Cz = 22,
                 output_file = "out_Backwater",
                 rho_w = 1000,
                 rho_s = 2650,
                 g = 9.81):

        self._friction = friction
        self._output_file = output_file

        self._qw = float(qw)
        self._initdepth = float(initdepth)
        self._initS = float(S)

        self._alpha = float(alpha)
        self._n = float(n)
        self._D50 = float(D50) / 1000 # M in old files
        self._D90 = float(D90) / 1000

        self._Cz = float(Cz) 

        self._g = float(g)

        self._L = int(L)
        self._dx = int(dx)


        ##### DERIVATIVES ######

        self._R = (float(rho_s) - float(rho_w)) / float(rho_w)

        self._nodes = (self._L / self._dx) + 1

        self._U = np.zeros((self._nodes,), dtype = np.float)
        self._taub = np.zeros_like(self._U)
        self._eta = np.zeros_like(self._U)
        self._H = np.zeros_like(self._U)
        self._x = np.zeros_like(self._U)
        self._S = np.zeros_like(self._U)
        
        self.set_surface()

        
    @property
    def bed_elevation(self):    
        return self._eta
        
    @bed_elevation.setter
    def bed_elevation(self, new_eta):
        
        msg = "The length of the bed elevation array " \
        "should match the length of the distance array (" + str(self._nodes) + " nodes)"  
        
        assert len(new_eta) == self._nodes, msg
        self._eta = new_eta
    
    @property
    def dam_height(self):
        """
        Returns the water surface elevation at the downstream boundary
        """
        return self._initdepth
    
    @dam_height.setter
    def dam_height(self, new_initdepth):
        self._initdepth = new_initdepth
    
    
    @property
    def discharge(self):
        """
        Returns the width-average water discharge at the upstream end
        """
        return self._qw
    
    @discharge.setter
    def discharge(self, new_qw):
        self._qw = new_qw
        
    
    @property
    def water_surface_elevation(self):
        """
        Returns the calculated water surface elevation
        """
        return self._ksi
    
    @property
    def water_depth(self):
        """
        Returns the calculated water depth
        """
        return self._H
    
    @property
    def flow_velocity(self):
        """
        Returns the calculated flow velocity
        """
        return self._U
    
    @property
    def bed_shear_stress(self):
        """
        Returns the calculated bed shear stress
        """
        return self._taub
    
    @property
    def distance(self):
        """
        Returns the x (distance) array
        """
        return self._x
    
    @property
    def friction_method(self):
        """
        Returns the friction method [Chezy, Manning]
        """
        return self._friction
    
    @friction_method.setter
    def friction_method(self, new_frict):
        """
        Set the friction method [Chezy, Manning]
        """
        
        friction_methods = ['Chezy', 'Manning']
        
        msg = "Friction method not recognized."
        
        assert new_frict in friction_methods, msg
        
        self._friction = new_frict
    



    ################################################
    ################# Functions ####################
    ################################################ 

    ################# Chezy ##################

    def Chezy_formulation(self):

        assert all([self._Cz]), \
            "One or more variables required for the Chezy formulation are missing"

        self.Hf = (self._qw**2 / (self._Cz**2 * self._g * self._initS))**(1./3)
        self.Uf = self._qw / self.Hf
        self.Fr = self.Uf  / sqrt(self._g * self.Hf)
        self.taubf = 1000 * (self.Uf**2 / self._Cz**2)

    def Chezy_Cf(self, Hm_):
        return 1 / self._Cz**2


    ########### Manning-Strickler ############

    def Manning_formulation(self):

        assert all([self._n, self._D90, self._D50, self._alpha]), \
            "One or more variables required for the Manning formulation are missing"

        ks = self._n * self._D90
        self.Hf = ((ks**(1./3) * self._qw**2) / \
                   (self._alpha**2 * self._g * self._initS))**(1./3)
        self.Uf = self._qw / self.Hf
        self.Fr = self.Uf  / sqrt(self._g * self.Hf)
        self.taubf = ((ks**(1./3) * self._qw**2) / \
                      (self._alpha**2 * self._g))**(1./3) * self._initS**0.7 / \
                      (self._R * self._D50)

            
    def Manning_Cf(self, Hm_):
        return ((self._n * self._D90) / Hm_)**(1./3) / self._alpha**2


    ########### Shared functions #############

    def set_surface(self):
        """
        Create an initial linear bed slope, in case the user doesn't define the surface
        """
        
        self._x[0] = self._L
        self._eta[0] = 0
        
        for i in range(1,self._nodes):
            self._x[i] = self._x[i-1] - self._dx
            self._eta[i] = self._eta[i-1] + (self._initS * self._dx)
    
    def get_critical_values(self):
        
        if self.friction_ == "Manning":
            self.Manning_formulation()

        if self.friction_ == "Chezy":
            self.Chezy_formulation()

        self.Hcrit = (self._qw**2 / self._g) ** (1./3)
        self.Ucrit = self._qw / self._initdepth

        
    def get_slope_array(self):
        self._S = -1*tools.central_difference(self._eta, self._dx)
        
        
    def check_direction(self):
        """
        Make sure that the _eta array increases with index
        (so the "dam" is at _x[0] and integration is upstream)
        """
        
        if self._eta[0] > self._eta[-1]:
            print 'flipping arrays'
            self._eta = np.flipud(self._eta)
            self._x = np.flipud(self._x)
            
    ###################################################
    
    
    def calculate_FnH(self, Cf_, H_, S_):

        FnH_ = ( S_ - ( Cf_ * self._qw**2 / ( self._g * H_**3 ))) / \
                ( 1 - ( self._qw**2 / ( self._g * H_**3 )))

        return FnH_
    
    
    def update_vals_Fn(self, Cf_fun, i):

        Cf = Cf_fun(self._H[i])
        FnH = self.calculate_FnH(Cf, self._H[i], self._S[i])

        Hp = self._H[i] - FnH * self._dx

        Cf = Cf_fun(Hp)
        FnHp = self.calculate_FnH(Cf, Hp, self._S[i])

        return FnH, FnHp, Cf

    
    def update_vals(self, i, Cf_):

        self._U[i] = self._qw / self._H[i]
        self._taub[i] = 1000 * Cf_ * self._U[i]**2


    def update_surface(self, i, FnH, FnHp):

        self._H[i] = self._H[i-1] - (0.5 * (FnH + FnHp) * self._dx)


    def Backwater_Calculator(self, Cf_fun):

        self._H[0] = self._initdepth
        FnH, FnHp, Cf = self.update_vals_Fn(Cf_fun, 0)
        self.update_vals(0, Cf)

        for i in range(1, self._nodes):

            self.update_surface(i, FnH, FnHp)
            FnH, FnHp, Cf = self.update_vals_Fn(Cf_fun, i)
            self.update_vals(i, Cf)
            
            
    ##########################################
    ################# Run ####################
    ##########################################
    
    def run(self):
        
        self.check_direction()
        self.get_slope_array()

        if self._friction == "Manning":
            self.Backwater_Calculator(self.Manning_Cf)

        if self._friction == "Chezy":
            self.Backwater_Calculator(self.Chezy_Cf)
            
        self._ksi = self._eta + self._H
            
     