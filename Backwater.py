import numpy as np
from math import sqrt
import json

class Backwater(object):
    
    def __init__(self,
                 friction = "Chezy",
                 qw = 5.7,
                 initdepth = 30.,
                 S = 0.00025,
                 total_x = 123000,
                 dx = 1000,
                 alpha = 8.1,
                 n = 2,
                 D90 = .91,
                 D50 = .27,
                 Cz = 22,
                 rho_w = 1000,
                 rho_s = 2650,
                 g = 9.81,
                 verbose = False,
                 save_output = True):
        

        self._verbose = verbose
        self.save_output = save_output
        
        self.friction = friction

        self.qw = float(qw)
        self.initdepth = float(initdepth)
        self.S = float(S)

        self.alpha = float(alpha)
        self.n = float(n)
        self.D50 = float(D50) / 1000 # M in old files
        self.D90 = float(D90) / 1000

        self.Cz = float(Cz) 

        self.g = float(g)
        self.rho_s = float(rho_s)
        self.rho_w = float(rho_w)

        self.total_x = int(total_x)
        self.dx = int(dx)


        ##### DERIVATIVES ######

        self.R = (self.rho_s - self.rho_w) / self.rho_w

        self.nodes = (self.total_x / self.dx) + 1

        self.U = np.zeros((self.nodes,), dtype = np.float)
        self.taub = np.zeros_like(self.U)
        self.eta = np.zeros_like(self.U)
        self.ksi = np.zeros_like(self.U)
        self.H = np.zeros_like(self.U)
        self.x = np.zeros_like(self.U)
        
        # Outputs
        self.Hf = None
        self.Uf = None
        self.Fr = None
        self.taubf = None
        
    @property
    def depth(self):
        return self.Hf

    @property
    def flow_velocity(self):
        return self.Uf
    
    @property
    def Froude_number(self):
        return self.Fr
    
    @property
    def shear_stress(self):
        return self.taubf
    
    
    ################################################
    ################# Functions ####################
    ################################################ 

    ################# Chezy ##################

    def Chezy_formulation(self):

        assert all([self.Cz]), \
            "One or more variables required for the Chezy formulation are missing"

        self.Hf = ( self.qw**2 / ( self.Cz**2 * self.g * self.S ) )**(1./3)
        
        self.Uf = self.qw / self.Hf
        
        self.Fr = self.Uf / sqrt( self.g * self.Hf )
        
        self.taubf = 1000. * ( self.Uf**2 / self.Cz**2 )

        
        
    def Chezy_Cf(self, Hm_):

        return 1 / self.Cz**2


    ########### Manning-Strickler ############

    def Manning_formulation(self):

        assert all([self.n, self.D90, self.D50, self.alpha]), \
            "One or more variables required for the Manning formulation are missing"

        ks = self.n * self.D90
        
        self.Hf = (( ks**(1./3) * self.qw**2 ) / \
                   ( self.alpha**2 * self.g * self.S ))**(1./3)
        
        self.Uf = self.qw / self.Hf
        
        self.Fr = self.Uf / sqrt( self.g * self.Hf )
        
        
        self.taubf = (( ks**(1./3) * self.qw**2 ) / \
                      ( self.alpha**2 * self.g ))**(1./3) * self.S**0.7 / \
                      ( self.R * self.D50)

            
            
    def Manning_Cf(self, Hm):

        return (( self.n * self.D90 ) / Hm )**(1./3) / self.alpha**2


    ########### Shared functions #############

    def get_critical_values(self):

        self.Hcrit = (self.qw**2 / self.g) ** (1./3)
        self.Ucrit = self.qw / self.initdepth
    
    
    def calculate_FnH(self, Cf, H):

        FnH = ( self.S - ( Cf * self.qw**2 / ( self.g * H**3 ))) / \
                ( 1 - ( self.qw**2 / ( self.g * H**3 )))

        return FnH
    
    
    def update_vals_Fn(self, Cf_fun, i):

        Cf = Cf_fun(self.H[i])
        FnH = self.calculate_FnH(Cf, self.H[i])

        Hp = self.H[i] - FnH * self.dx

        Cf = Cf_fun(Hp)
        FnHp = self.calculate_FnH(Cf, Hp)

        return FnH, FnHp, Cf

    
    def update_vals(self, i, Cf):

        self.U[i] = self.qw / self.H[i]
        self.taub[i] = 1000. * Cf * self.U[i]**2
        self.ksi[i] = self.H[i] + self.eta[i]


    def update_surface(self, i, FnH, FnHp):

#         self.x[i] = self.x[i-1] - self.dx_
        self.H[i] = self.H[i-1] - (0.5 * (FnH + FnHp) * self.dx)
#         self.eta[i] = self.eta[i-1] + (self.S_ * self.dx_)


    def Backwater_Calculator(self, Cf_fun):

#         self.x[0] = self.total_x_
        self.H[0] = self.initdepth
#         self.eta[0] = 0
        FnH, FnHp, Cf = self.update_vals_Fn(Cf_fun, 0)
        self.update_vals(0, Cf)

        for i in range(1, self.nodes):

            self.update_surface(i, FnH, FnHp)
            FnH, FnHp, Cf = self.update_vals_Fn(Cf_fun, i)
            self.update_vals(i, Cf)

    def set_surface(self):
        
        self.x[0] = self.total_x
        self.eta[0] = 0
        
        for i in range(1,self.nodes):
            self.x[i] = self.x[i-1] - self.dx
            self.eta[i] = self.eta[i-1] + (self.S * self.dx)
            
            
    ##########################################
    ################# Run ####################
    ##########################################   
    def run(self): 

        self.get_critical_values()
        
        self.set_surface()

        if self.friction == "Manning":
            self.Manning_formulation()
            self.Backwater_Calculator(self.Manning_Cf)


        if self.friction == "Chezy":
            self.Chezy_formulation()
            self.Backwater_Calculator(self.Chezy_Cf)
            
     
    
    ##########################################
    ############## Finalize ##################
    ##########################################   
    def finalize(self):
        
        if self.save_output:
            
            output_dict = {
                'Depth' : self.Hf,
                'Flow_velocity' : self.Uf,
                'Froude_number' : self.Fr,
                'Shear_stress' : self.taubf
            }
            
            with open('output/Backwater_steadystate.json', 'w') as f:
                json.dump(output_dict, f, indent=4)
                


            fields = ['Position', 'Depth', 'FlowVelocity',
                      'ShearStress', 'BedElevation', 'WaterSurfaceElev']

            header = ', '.join(fields)

            dtype = [(i, float) for i in fields]
            self.data = np.empty(len(self.x), dtype = dtype)

            self.data['Position'] = self.x
            self.data['Depth'] = self.H
            self.data['FlowVelocity'] = self.U
            self.data['ShearStress'] = self.taub
            self.data['BedElevation'] = self.eta
            self.data['WaterSurfaceElev'] = self.ksi

            np.savetxt('output/Backwater_profile.csv', self.data,
                       header = header,
                       delimiter = ",",
                       fmt = '%10.5f',
                       comments = '')

            if self._verbose:
                print 'Saved data to file'