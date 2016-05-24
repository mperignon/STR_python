import numpy as np
from math import sqrt, log
import shared_tools as tools

from collections import OrderedDict
import json

class UnsortableList(list):
    def sort(self, *args, **kwargs):
        pass
    
class UnsortableOrderedDict(OrderedDict):
    def items(self, *args, **kwargs):
        return UnsortableList(OrderedDict.items(self, *args, **kwargs))


class DeltaBW(object):
    """
    DeltaBW (1DDeltaBW)

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)
    

    MC Perignon
    Nov 2015

    ----------------------------------------------

    Computes the evolution of the long profile of a river ending in a 1D
    migrating delta using a backwater formulation

    """

    
    def __init__(self,
                friction = "Chezy",
                flood_discharge = 6.,
                flood_intermittency = 1.,
                bedload_discharge_input = 0.001,
                grain_size = 0.5,
                coeff_Chezy = 15.,
                exponent_load_relation = 2.5,
                critical_Shields_stress = 0.,
                foreset_elev__top__init = 3.,
                foreset_elev__base__init = 0.,
                bed_slope__init = 0.00025,
                basin_slope = 0.,
                domain_length = 10000.,
                domain_length_max = 500000.,
                foreset_slope = 0.2,
                bed_porosity = 0.4,
                coeff_Manning_ks = 0.3,
                coeff_bedload = 7.2,
                coeff_Manning = 8.1,
                dx = 500.,
                dt__days = 0.182625,
                num_iterations = 30000,
                print_step = 5000,
                density_sediment = 2650.,
                density_water = 1000.,
                basin_water_surface_elevation = 8.5,
                verbose = False,
                save_output = True):

        
        self.verbose = verbose
        self.save_output = save_output
        
        day_to_sec = 60. * 60. * 24.
        self.sec_to_year = 1/(365.25 * 24 * 60 * 60)
        self.g = 9.81
        self.R = (float(density_sediment) - float(density_water)) / float(density_water)

        self.friction = friction
        self.flood_discharge = float(flood_discharge)
        self.flood_intermittency = float(flood_intermittency)
        self.bedload_discharge_input = float(bedload_discharge_input)
        self.grain_size = grain_size / 1000.
        self.coeff_Chezy = float(coeff_Chezy)
        self.exponent_load_relation = float(exponent_load_relation)
        self.critical_Shields_stress = float(critical_Shields_stress)
        self.basin_slope = float(basin_slope)
        self.domain_length = float(domain_length)
        self.domain_length_max = float(domain_length_max)
        self.foreset_slope = float(foreset_slope)
        self.bed_porosity = float(bed_porosity)
        self.coeff_Manning_ks = float(coeff_Manning_ks) / 1000.
        self.coeff_bedload = float(coeff_bedload)
        self.coeff_Manning = float(coeff_Manning)
        self.dx = float(dx)
        self.dt = float(dt__days) * day_to_sec
        self.num_iterations = int(num_iterations)
        self.print_step = int(print_step)
        
        self.basin_water_surface_elevation = float(basin_water_surface_elevation)
        
        


        num_nodes = int(self.domain_length / self.dx) + 2

        self.dxbar = 1. / (num_nodes - 2)
        self.xfbar = self.dxbar * np.arange(num_nodes)
        self.foreset_elev__top = foreset_elev__top__init
        self.foreset_elev__base = foreset_elev__base__init

        self.eta_upstream = (self.domain_length * bed_slope__init +
                            self.foreset_elev__top)

        self.bed_length = (self.domain_length +
                    (self.foreset_elev__top - self.foreset_elev__base) / self.foreset_slope)



        self.x = np.zeros((num_nodes,))
        self.x = self.domain_length * self.xfbar
        self.x[-1] = self.bed_length

        self.eta = np.zeros_like(self.x)
        self.eta = (self.domain_length * bed_slope__init -
                    self.domain_length * bed_slope__init * self.xfbar +
                    self.foreset_elev__top)
        self.eta[-1] = self.foreset_elev__base

        self.bedload_discharge__div = np.zeros_like(self.x)
        self.eta__div = np.zeros_like(self.x)
        
        self.Hfluv = np.zeros_like(self.x)
        self.etafluv = self.eta.copy()
        
        self.sed_volume__init = self.calculate_mass_balance()
        
        # Output
        num_prints = int(num_iterations / self.print_step)
        
        self.fields = [str(i*self.dt*self.sec_to_year) + '_years' for i in range(0,num_iterations+1,int(print_step))]

        dtype = [(i, float) for i in self.fields]
        self.print_x = np.empty((num_nodes,), dtype = dtype)
        self.print_eta = np.empty((num_nodes,), dtype = dtype)
        self.print_qb = np.empty((num_nodes,), dtype = dtype)
        self.print_H = np.empty((num_nodes,), dtype = dtype)
        self.print_tau = np.empty((num_nodes,), dtype = dtype)
        
        self.output_dict = UnsortableOrderedDict()
        


    def run(self):
        
        self.compute_normal_flow()

        for j in range(self.num_iterations+1):

            time = self.dt * j * self.sec_to_year

            self.get_water_depth()

            self.get_bedload_discharge()

            self.solve_Exner_equation()


            if (j) % self.print_step == 0:
                
                self.record_output(j)
                
                if self.verbose:
                    print 'time:', time 

                
    def record_output(self, j):
        
        ind = self.fields[int(j / self.print_step)]
        
        self.print_x[ind] = self.x
        self.print_eta[ind] = self.eta
        self.print_qb[ind] = self.bedload_discharge
        self.print_H[ind] = self.Hfluv
        self.print_tau[ind] = self.shear_stress
        
        self.output_dict[ind] = {
            
            'domain_length__sss' : self.domain_length,
            'bed_length__sbb' : self.bed_length,
            'bed_elev_upstream__etaup' : self.eta_upstream,
            'bed_elev_foreset_top__etatop' : self.eta[-2],
            'bed_elev_foreset_base__etabottom' : self.eta[-1]
            
            }
        
        
  

    def compute_normal_flow(self):
        
        tsn = ((self.bedload_discharge_input / (self.coeff_bedload * self.grain_size *
            sqrt(self.R * self.g * self.grain_size)))**(1. / self.exponent_load_relation) +
               self.critical_Shields_stress)
       
       
        if self.friction == 'Chezy':

            Cf = 1. / (self.coeff_Chezy**2)
            slope__normal = (self.coeff_Chezy  * (self.R * self.grain_size * tsn)**1.5 *
                 sqrt(self.g) / self.flood_discharge)

        if self.friction == 'Manning':
            
            slope__normal = ((tsn * self.R * self.grain_size)**(10./7) *
                 (self.coeff_Manning**2 * self.g /
                  (self.coeff_Manning_ks**(1./3) * self.flood_discharge**2))**(3./7))


        self.depth__normal = self.R * self.grain_size * tsn / slope__normal

        Froude__normal = (self.flood_discharge / self.depth__normal) / sqrt(self.g * self.depth__normal)

        assert Froude__normal < 1, 'The Froude number for normal flow is greater than 1'

        Hi = self.basin_water_surface_elevation - self.foreset_elev__top
        Hcrit = (self.flood_discharge**2 / self.g)**(1./3)

        assert (Hi > Hcrit), 'Hi must be greater than Hcrit for backwater calculations'

        
        

    def get_water_depth(self):
        
        self.Hfluv[-1] = self.basin_water_surface_elevation - self.etafluv[-2]

        if self.friction == 'Chezy':
            
            Cf = 1. / (self.coeff_Chezy ** 2)
            
            for i in range(len(self.x)-1):

                fr2p = self.flood_discharge**2 / (self.g * self.Hfluv[-1-i]**3)
                
                fnp = ((self.etafluv[-2-i] - self.etafluv[-1-i] -
                        Cf * fr2p * self.domain_length * self.dxbar) / (1. - fr2p))
                       
                Hpred = self.Hfluv[-1-i] - fnp
                
                fr2 = self.flood_discharge**2 / (self.g * Hpred ** 3)
                fn = ((self.etafluv[-2-i] - self.etafluv[-1-i] -
                        Cf * fr2 * self.domain_length * self.dxbar) / (1. - fr2))

                self.Hfluv[-2-i] = self.Hfluv[-1-i] - 0.5 * (fn + fnp)

                if self.Hfluv[-2-i] < self.depth__normal:
                    self.Hfluv[-2-i] = self.depth__normal

        if self.friction == 'Manning':
            
            for i in range(len(x)-1):
                
                Cf = (1.0 / self.coeff_Manning**2) * (self.coeff_Manning_ks / self.Hfluv[-1-i])**(1./3)

                fr2p = self.flood_discharge**2 / (self.g * self.Hfluv[-1-i]**3)
                
                fnp = ((self.etafluv[-2-i] - self.etafluv[-1-i] -
                        Cf * fr2p * self.domain_length * dxbar) / (1. - fr2p))
                       
                Hpred = self.Hfluv[-1-i] - fnp
                
                fr2 = self.flood_discharge**2 / (self.g * Hpred ** 3)
                fn = ((self.etafluv[-2-i] - self.etafluv[-1-i] -
                        Cf * fr2 * self.domain_length * dxbar) / (1. - fr2))

                self.Hfluv[-2-i] = self.Hfluv[-1-i] - 0.5 * (fn + fnp)

                if self.Hfluv[-2-i] < self.depth__normal:
                    self.Hfluv[-2-i] = self.depth__normal



    def get_bedload_discharge(self):
        
        if self.friction == 'Chezy':
            
            Cf = 1. / (self.coeff_Chezy ** 2)
        
        if self.friction == 'Manning':
            
            Cf = (1.0 / self.coeff_Manning**2) * (self.coeff_Manning_ks / self.Hfluv)**(1./3)

            
        self.shear_stress = Cf * (self.flood_discharge / self.Hfluv)**2 / (self.R * self.g * self.grain_size)
        
        self.bedload_discharge = np.zeros_like(self.x)
        
        self.bedload_discharge[self.shear_stress > self.critical_Shields_stress] = (
            sqrt(self.R * self.g * self.grain_size) * self.grain_size * self.coeff_bedload *
            (self.shear_stress[self.shear_stress > self.critical_Shields_stress] -
             self.critical_Shields_stress)**self.exponent_load_relation)



    def solve_Exner_equation(self):
        
        domain_length__change = (self.flood_intermittency * self.bedload_discharge[-2] /
                (self.foreset_slope * (1 - self.bed_porosity) *
                (self.bed_length - self.domain_length)))
        
        etatopold = self.foreset_elev__top


        self.bedload_discharge__div[0] = (self.bedload_discharge[1] - self.bedload_discharge_input) / self.dxbar

        self.eta__div[0] = (self.etafluv[1] - self.etafluv[0]) / self.dxbar


        for i in range(1, len(self.eta)-1):

            self.bedload_discharge__div[i] = (self.bedload_discharge[i] - self.bedload_discharge[i-1]) / self.dxbar

            self.eta__div[i] = (self.etafluv[i+1] - self.etafluv[i-1]) / (2 * self.dxbar)


        self.bedload_discharge__div[-1] = (self.bedload_discharge[-1] - self.bedload_discharge[-2]) / self.dxbar
        self.eta__div[-1] = (self.etafluv[-1] - self.etafluv[-2]) / self.dxbar

        self.etafluv += (self.dt * (- self.flood_intermittency * self.bedload_discharge__div /
                (self.domain_length * (1 - self.bed_porosity)) +
                self.xfbar * domain_length__change * self.eta__div / self.domain_length))

        self.eta = self.etafluv.copy()
        

        self.foreset_elev__top = self.etafluv[-2]
        
        etasdot = (self.foreset_elev__top - etatopold) / self.dt
        domain_length__change +=  - etasdot / self.foreset_slope

        self.domain_length += domain_length__change * self.dt
        

        self.x = self.domain_length * self.xfbar
        self.x[-1] = self.bed_length
        
        
        bed_length__change = ((self.foreset_slope * domain_length__change + etasdot) /
                              (self.foreset_slope - self.basin_slope))
        
        self.bed_length += bed_length__change * self.dt


        self.eta_upstream = self.eta[0]
        self.foreset_elev__base += -self.basin_slope * bed_length__change / self.dt

        if self.foreset_elev__base < 1e-6:
            self.foreset_elev__base = 0

        self.eta[-1] = self.foreset_elev__base

        assert self.bed_length < self.domain_length_max, 'The delta prograded past the max length of the basin'
        
    
    
    def calculate_mass_balance(self):

        abase = (self.foreset_elev__base + 
                 self.basin_slope * (self.domain_length + 
                (self.foreset_elev__top - self.foreset_elev__base) / self.foreset_slope))

        etabase = abase - self.basin_slope * self.x
        etad = self.eta - etabase

        sediment_volume = 0.5 * etad[0] * (self.x[1] - self.x[0])
        for i in range(len(etad)-1):
            sediment_volume += 0.5 * (etad[i] + etad[i+1]) * (self.x[i+1] - self.x[i])

        sediment_volume += 0.5 * etad[-1] * (self.x[-1] - self.x[-2])
        sediment_volume = sediment_volume * (1.0 - self.bed_porosity)
        
        return sediment_volume
    
    
    
    def finalize(self):
        
        self.sed_volume__final = self.calculate_mass_balance()
        self.sed_volume__feed = self.bedload_discharge_input * self.dt * self.num_iterations
        
        self.mass_balance_error = abs(
            (self.sed_volume__final - self.sed_volume__init - self.sed_volume__feed) / 
             self.sed_volume__final)
        
        
        
        if self.save_output:
        
            header = ', '.join(['x','eta','qb','H','tau'])

            for i in range(len(self.fields)):

                self.data = np.zeros((len(self.x),5))
                self.data[:,0] = self.print_x[self.fields[i]]
                self.data[:,1] = self.print_eta[self.fields[i]]
                self.data[:,2] = self.print_qb[self.fields[i]]
                self.data[:,3] = self.print_H[self.fields[i]]
                self.data[:,4] = self.print_tau[self.fields[i]]

                np.savetxt('output/DeltaBW_' + self.fields[i] + '.csv', self.data,
                               header = header,
                               delimiter = ",",
                               fmt = '%10.5f',
                               comments = '')

            self.output_dict['mass_balance'] = {

                'sed_volume_initial' : self.sed_volume__init,
                'sed_volume_final' : self.sed_volume__final,
                'sed_volume_total_feed' : self.sed_volume__feed,
                'mass_balance_error' : self.mass_balance_error

                }


            with open('output/DeltaBW_stats.json', 'w') as f:
                json.dump(self.output_dict, f, indent=4)
        
        
        
        