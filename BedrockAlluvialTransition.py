import numpy as np
from math import pow, sqrt

class BedrockAlluvialTransition(object):
    """
    BedrockAlluvialTransition

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Jan 2016

    ----------------------------------------------

    Calculates aggradation and degradation of a 1D profile with a migrating
    bedrock-alluvial transition at the upstream end. The bedrock basement channel
    is assumed to undergo no incision on the timescales at which the alluvial reach
    responds to change. The normal flow approximation is used and base level is
    maintained at eta == 0. The Engelund-hansen relation is used to compute sediment
    transport rates so the analyis is appropriate for sand-bed streams.

    """    
    def __init__(self,
                 flood_discharge = 6.,
                 flood_intermittency = 0.1,
                 sed_discharge = 0.0006,
                 grain_size = 0.5,
                 Chezy_coeff = 15.,
                 slope__bedrock = 0.00035,
                 slope__alluvium_init = 0.00009,
                 porosity = 0.4,
                 reach_length = 50000.,
                 dx = 2500.,
                 dt = 0.005,
                 total_time = 250,
                 print_dt = 50,
                 critical_shear_stress = 0.,
                 rho_w = 1000.,
                 rho_s = 2650.,
                 verbose = False,
                 save_output = True):
        
        self.g = 9.81
        self.ns = 2.5
        
        self.verbose = verbose
        self.save_output = save_output

        self.yr_to_sec = 31557600

        self.dt_sec = dt * self.yr_to_sec
        self.dt = dt
        self.total_time = total_time * self.yr_to_sec
        self.print_dt = print_dt * self.yr_to_sec

        self.flood_discharge = float(flood_discharge)
        self.flood_intermittency = float(flood_intermittency)
        self.sed_discharge = float(sed_discharge)
        self.Chezy_coeff = float(Chezy_coeff)
        self.slope__bedrock = float(slope__bedrock)
        self.slope__alluvium_init = float(slope__alluvium_init)
        self.porosity = float(porosity)
        self.reach_length = float(reach_length)
        self.dx = float(dx)
        self.tau_crit = float(critical_shear_stress)
        
        self.grain_size = float(grain_size) / 1000.

        self.R = (float(rho_s) - float(rho_w)) / float(rho_w)
        self.alpha = 0.05 * self.Chezy_coeff ** 2

        num_nodes = int(self.reach_length / self.dx) + 1

        # arrays

        self.xbar = np.arange(0,num_nodes) / float(num_nodes - 1)
        self.x = np.arange(0, self.reach_length + 1, self.dx)
        self.slope = np.ones_like(self.x) * self.slope__alluvium_init
        self.eta = (self.reach_length * self.slope__alluvium_init - 
                    self.reach_length * self.slope__alluvium_init * self.xbar)

        self.bedload_discharge = np.zeros_like(self.x)
        self.depth = np.zeros_like(self.x)
        self.shear_stress = np.zeros_like(self.x)


        self.x_transition = 0.
        self.dxbar = 1. / (num_nodes - 1)
        
        
        # Output
        num_prints = int(self.total_time / self.print_dt)
        
        max_time = int(self.total_time) / self.yr_to_sec
        d_time = int(self.print_dt) / self.yr_to_sec
        
        self.fields = [str(i) + '_years' for i in range(0,max_time+1, d_time)]

        dtype = [(i, float) for i in self.fields]
        self.print_x = np.empty((len(self.x),), dtype = dtype)
        self.print_eta = np.empty((len(self.x),), dtype = dtype)
        self.print_depth = np.empty((len(self.x),), dtype = dtype)
        self.print_slope = np.empty((len(self.x),), dtype = dtype)
        self.print_qb = np.empty((len(self.x),), dtype = dtype)

        
    @property
    def bed_elevation(self):
        return self.eta
    
    @property
    def water_surface_elevation(self):
        return self.eta + self.depth
    
    @property
    def water_depth(self):
        return self.depth
    
    @property
    def alluvial_slope(self):
        return self.slope
        
    @property
    def bed_shear_stress(self):
        return self.shear_stress
        
    @property
    def distance(self):
        return self.x
    
    @property
    def bedload_sed_discharge(self):
        return self.bedload_discharge
    
    @property
    def bedrock_alluvial_transition(self):
        return self.x_transition

    def run(self):
                    
        time = 0

        while time <= self.total_time:

            self.dx_ = self.dxbar * (self.reach_length - self.x_transition)
            
            self.find_bedload()
            self.find_new_eta()
            
            if time % self.print_dt == 0:
                
                self.record_output(time)

            time = time + self.dt_sec
            
            
                            
                
                
                
    def record_output(self, t):
        
        ind = self.fields[int(t / self.print_dt)]
        
        self.print_x[ind] = self.x
        self.print_eta[ind] = self.eta
        self.print_depth[ind] = self.depth
        self.print_slope[ind] = self.slope
        self.print_qb[ind] = self.bedload_discharge
        
        
        
    def finalize(self):
        
        if self.save_output:
        
            header = ', '.join(['x','eta','depth','slope','bedload_discharge'])

            for i in range(len(self.fields)):

                self.data = np.zeros((len(self.x),5))
                self.data[:,0] = self.print_x[self.fields[i]]
                self.data[:,1] = self.print_eta[self.fields[i]]
                self.data[:,2] = self.print_depth[self.fields[i]]
                self.data[:,3] = self.print_slope[self.fields[i]]
                self.data[:,4] = self.print_qb[self.fields[i]]

                np.savetxt('output/BedrockAlluvialTransition_' + self.fields[i] + '.csv', self.data,
                               header = header,
                               delimiter = ",",
                               fmt = '%10.8f',
                               comments = '')
            
            
    def find_bedload(self):
        
        St = np.zeros_like(self.x)
        self.slope[0] = (self.eta[0] - self.eta[1]) / self.dx_

        for i in range(1, len(self.slope)-1):
            self.slope[i] = (self.eta[i-1] - self.eta[i+1]) / (2 * self.dx_)
            St[i] = (self.eta[i-1] - self.eta[i]) / self.dx_

        St[0] = St[1]
        St[-1] = (self.eta[-2] - self.eta[-1]) / self.dx_


        self.depth = (self.flood_discharge**2 /
                      (self.g * St * self.Chezy_coeff**2))**(1./3)
        self.depth[0] = self.depth[1]

        self.shear_stress = self.depth * St / (self.R * self.grain_size)

        tau_diff = self.shear_stress - self.tau_crit
        tau_diff[tau_diff < 0] = 0
        tau_pow = np.array(map(lambda x: pow(x, self.ns), tau_diff))

        self.bedload_discharge = (self.alpha * tau_pow *
                                  sqrt(self.R * self.g * self.grain_size) *
                                  self.grain_size)

        self.bedload_discharge[0] = self.sed_discharge
        
        
    def find_new_eta(self):
        
        eta_upstream_old = self.eta[0]
        self.eta[0] += ((self.slope__bedrock / (self.slope__bedrock - self.slope[0])) *
                        (self.bedload_discharge[0] - self.bedload_discharge[1]) *
                         self.dt_sec / self.dx_)

        x_transition_dot = ((eta_upstream_old - self.eta[0]) /
                           (self.slope__bedrock * self.dt_sec))

        for i in range(1, len(self.x)-1):

            eta_sum = (-self.dt_sec * self.slope[i] * x_transition_dot * (1 - self.xbar[i]) +
                      self.flood_intermittency * self.dt_sec *
                      (self.bedload_discharge[i] - self.bedload_discharge[i+1]) /
                      (self.dx_ * (1 - self.porosity)))

            self.eta[i] += eta_sum

        self.x_transition += x_transition_dot * self.dt_sec

        self.x = self.x_transition + self.xbar * (self.reach_length - self.x_transition)