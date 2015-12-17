import numpy as np
import json

import GSD_calculator as gsd

from math import sqrt, log, exp
from shared_tools import po, oo, so
import shared_tools as tools

import Acronym

default_GSD = [(32,100),
               (16,50),
               (8,0)]


default_flow = [(80,0),
                (40, 20),
                (30, 60),
                (10, 100)]

class Acronym_D(Acronym.Acronym):

    def __init__(self,
                 GSD = default_GSD,
                 flow_curve = default_flow,
                 channel_slope = 0.005,
                 channel_width = 12.,
                 mannings_n = 2.,
                 rho_s = 2650.,
                 rho_w = 1000.,
                 alpha_r = 8.1,
                 verbose = False,
                 save_output = True):

        self.verbose = verbose
        self.save_output = save_output
        
        self.GSD = GSD
        self.flow_edges = flow_curve
        self.channel_width = float(channel_width)
        self.channel_slope = float(channel_slope)
        self.mannings_n = float(mannings_n)
        
        self.R = (float(rho_s) - float(rho_w))/ float(rho_w)
        self.g = 9.81
        
        self.alpha_r = float(alpha_r)

        self._percent_finer = None
        self._characteristic_size = None
        self.tau_star_surface = None
        self._surface_statistics = None
        self._bedload_statistics = None
    
    
    @property
    def Shields_stress(self):
        
        return self.bedload_averages['Shields']
    
    @property
    def averages(self):
        
        return self.bedload_averages
    
    @property
    def discharge_distribution(self):
        
        return self.flow_edges
    
    @discharge_distribution.setter
    def discharge_distribution(self, new_qs):
        
        self.flow_edges = new_qs
        
        

    def run(self):

        self.GSD_edges, self.GSD_surface = self.process_surface_GSD(self.GSD)

        self.pf_surface = self.GSD_edges['fractions']*100
        
        
        if len(self.GSD_surface) > 1:
            
            surface_vals = self.grain_statistics(self.GSD_edges,
                                                 self.GSD_surface,
                                                 self.pf_surface)
            
            self.c_sizes_surface = surface_vals[0]
            self.D_g_surface = surface_vals[1]
            self.sigma_g_surface = surface_vals[2]
            self.sigma_surface = surface_vals[3]
            
            
            
        self.flow_curve = self.process_flow_curve(self.flow_edges)
            
        for i in range(len(self.flow_curve)):
            self.calculate_one_flow_value(i)
        
        self.bedload_averages = self.get_bedload_averages(self.flow_curve)
        
        

        if len(self.GSD_surface) > 1:
            
            self.pf_bedload = self.bedload_averages['percent_finer'][0]
            
            bedload_vals = self.grain_statistics(self.GSD_edges,
                                                 self.GSD_surface,
                                                 self.pf_bedload)
            
            self.c_sizes_bedload = bedload_vals[0]
            self.D_g_bedload = bedload_vals[1]
            self.sigma_g_bedload = bedload_vals[2]
            self.sigma_bedload = bedload_vals[3]
        


    def calculate_one_flow_value(self, i):
        
        Q = self.flow_curve['Q'][i]
        pQ = self.flow_curve['pQ'][i]


        ustar, flow_depth = self.get_flow_characteristics(Q)

        self.flow_curve['ustar'][i] = ustar
        self.flow_curve['H'][i] = flow_depth


        calc_bedload_vals = self.calculate_bedload(self.GSD_surface,
                                               self.D_g_surface,
                                               self.sigma_surface,
                                               ustar)

        self.flow_curve['Shields'][i] = calc_bedload_vals[0]
        self.flow_curve['q_bedload'][i] = calc_bedload_vals[3] 


        GSD_bedload = calc_bedload_vals[1]
        pf_bedload = self.get_percent_finer(GSD_bedload)

        self.flow_curve['percent_finer'][i] = pf_bedload
        
        
        
        
    def get_bedload_averages(self, flow_curve):
        
        type_pf = str(len(self.GSD_edges)) + 'float'
        dtype = [('Q', float), ('H', float), ('ustar', float),
                 ('Shields', float), ('q_bedload', float), ('percent_finer', type_pf)]
        averages = np.zeros((1), dtype=dtype)
        
        
        averages['Q'] = sum(flow_curve['Q'] * flow_curve['pQ'])
        averages['H'] = sum(flow_curve['H'] * flow_curve['pQ'])
        averages['ustar'] = sum(flow_curve['ustar'] * flow_curve['pQ'])
        averages['Shields'] = sum(flow_curve['Shields'] * flow_curve['pQ'])
        averages['q_bedload'] = sum(flow_curve['q_bedload'] * flow_curve['pQ'])
        averages['percent_finer'] = self.flow_curve['pQ'].dot(self.flow_curve['percent_finer'])
            
        return averages
        

        
    
            
    def process_flow_curve(self, flow_edges):
        
        # get the middle value on the interval
        
        type_pf = str(len(self.GSD_edges)) + 'float'

        dtype = [('Q', float), ('pQ', float), ('H', float), ('ustar', float),
                 ('Shields', float), ('q_bedload', float), ('percent_finer', type_pf)]
        flow_curve = np.zeros((max((len(self.flow_edges)-1, 1))), dtype=dtype)
        
        if len(flow_edges) == 1:
            
            flow_curve['Q'] = float(flow_edges[0][0])
            flow_curve['pQ'] = float(flow_edges[0][1]) / 100
            
        else:
            
            # get midpoints            
            for i in range(len(flow_edges)-1):
                
                flow_curve['Q'][i] = (float(flow_edges[i][0]) + float(flow_edges[i+1][0])) / 2
                flow_curve['pQ'][i] = (float(flow_edges[i+1][1]) - float(flow_edges[i][1])) / 100       

        return flow_curve        

    
    
    
    
    def get_flow_characteristics(self, Q):
        
        Dx90s = self.c_sizes_surface[-1]
        
        ks = self.mannings_n * (Dx90s / 1000)
        Ho = ((ks**(1./3) * Q**2) / 
              (self.alpha_r**2 * self.g * self.channel_width**2 * self.channel_slope))**(1./3)
        
        max_possible_iterations = 100
        error_threshold = 0.001
        H = Ho
        i = 0
        
        while i < 100:
            
            i += 1
            
            Hnew = Ho * (1 - (H / self.channel_width))**(-1./3)
            err = abs(2 * (Hnew - H)/(Hnew + H))
            H = Hnew
            
            if err <= error_threshold:
                 i = 1000 # break
                
        assert i != 100, 'Depth calculation failed'
        
        ustar = sqrt(self.g * H * self.channel_slope * (1 - (H / self.channel_width)))
        
        return ustar, H
    
    
    
    def finalize(self):
        
        if self.save_output:

            surface_dict = {

                'Geometric_mean' : self.surface_statistics['Dg'],
                'Standard_Deviation' : self.surface_statistics['sigma'],
                'Geometric_Standard_Deviation' : self.surface_statistics['sigma_g'],
                'D_30' : self.characteristic_size['surface'][0],
                'D_50' : self.characteristic_size['surface'][1],
                'D_70' : self.characteristic_size['surface'][2],
                'D_90' : self.characteristic_size['surface'][3]

                }

            bedload_dict = {

                'Geometric_mean' : self.bedload_statistics['Dg'],
                'Standard_Deviation' : self.bedload_statistics['sigma'],
                'Geometric_Standard_Deviation' : self.bedload_statistics['sigma_g'],
                'D_30' : self.characteristic_size['bedload'][0],
                'D_50' : self.characteristic_size['bedload'][1],
                'D_70' : self.characteristic_size['bedload'][2],
                'D_90' : self.characteristic_size['bedload'][3]

                }
            
            flow_dict = [dict(zip(self.averages.dtype.names, r)) for r in self.averages][0]
            p_finer_dict = [dict(zip(list(self.percent_finer['diameter']), flow_dict['percent_finer']))]
            flow_dict['percent_finer'] = p_finer_dict


            stats_output = {

                'Surface' : surface_dict,
                'Bedload' : bedload_dict,
                'Averages': flow_dict

                }

            with open('output/AcronymD_stats.json', 'w') as f:
                    json.dump(stats_output, f, indent=4)
        
        
            fields = ['GrainSize', 'PercentFiner_Surface', 'PercentFiner_Bedload']

            header = ', '.join(fields)

            np.savetxt('output/AcronymD_PercentFiner.csv', self.percent_finer,
                       header = header,
                       delimiter = ",",
                       fmt = '%10.5f',
                       comments = '')
        