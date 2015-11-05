import numpy as np
from math import sqrt, log
import matplotlib.pyplot as plt

import shared_tools as tools

class SteadyStateAg(object):
    
    def __init__(self,
                 
                 time = 10000.,
                 channel_length = 400000.,
                 floodplain_width = 8000.,
                 channel_sinuosity = 1.5,
                 bed_porosity = 0.4,
                 flood_intermittence = 0.1,
                 bankfull_discharge = 5000.,
                 annual_bedload_discharge = 12.,
                 grain_size = 0.25,
                 fraction_washload_deposition = 2.5,
                 sea_level_elevation = 0.,
                 sea_level_change = 0.010,
                 Shields_number = 1.86,
                 coeff_Engelund = 0.05,
                 coeff_Chezy = 25.,
                 density_water = 1000.,
                 density_sediment = 2650.,
                 gravity = 9.81):

        self.yr_to_sec = 60 * 60 * 24 * 365.25
        self.g = float(gravity)
        
        
        
        self.channel_length = float(channel_length) # m
        self.floodplain_width = float(floodplain_width) # m
        self.channel_sinuosity = float(channel_sinuosity)
        self.porosity = float(bed_porosity)
        
        self.flood_intermittence = float(flood_intermittence)
        self._bankfull_discharge = float(bankfull_discharge) # m3/s
        
        self.annual_bedload_discharge = (
            float(annual_bedload_discharge) *  1e6 / self.yr_to_sec) # Mt/yr to m3/s
        self.grain_size = float(grain_size) / 1000. # mm to m
        self.fraction_washload_deposition = float(fraction_washload_deposition)
        
        self.sea_level_elevation = float(sea_level_elevation)
        self.sea_level_change = float(sea_level_change) # m/yr
        
        self.Shields_number = float(Shields_number)
        self.coeff_Engelund = float(coeff_Engelund)
        self.coeff_Chezy = float(coeff_Chezy)
        
        self.time = float(time)
        

        self.R = (float(density_sediment) - float(density_water)) / float(density_water)

        # set up initial bed
        self.distance__nd = np.array(np.arange(0,1.05,0.05))
        self.distance = self.channel_length * self.distance__nd
        
        

        
    def run(self):

        self.Cf = tools.Chezy_Cf(self.coeff_Chezy)
    
        self._bedload_discharge = self.annual_bedload_discharge / ((self.R + 1.) * self.flood_intermittence)
        self._bedload_concentration = self._bedload_discharge / (self._bedload_discharge + self._bankfull_discharge)
    
    
        
        self.beta = (
            (1. - self.porosity) * self.floodplain_width * 
            (self.sea_level_change / self.yr_to_sec) *  self.channel_length /
            (self.flood_intermittence * self.channel_sinuosity * self._bedload_discharge *
            (1. + self.fraction_washload_deposition)))

        assert self.beta < 1, (
            "Adjust L such that the beta parameter is below 1. "
            "If beta is greater than 1, the sediment does not reach the delta.")

        
        self._bed_slope = (
            self.R * sqrt(self.Cf) * self._bedload_discharge /
            (self.coeff_Engelund * self.Shields_number * self._bankfull_discharge))
        
        
        self.bed_elevation_change__nd = (
            self._bed_slope * ((1 - 0.5 * self.beta) - self.distance__nd
                                 + 0.5 * self.beta * self.distance__nd**2))
        self.bed_elevation_change__nd[-1] = 0.
        
        self.bed_elevation_change = self.channel_length * self.bed_elevation_change__nd
        
        
        
        
    def get_bed_elevation(self, time):
        
        z = (self.sea_level_elevation
                        + self.sea_level_change * float(time)
                        + self.bed_elevation_change)
        
        return z
    
        
    
    @property
    def bed_slope(self):
        
        S = self._bed_slope * (1 - self.beta * self.distance__nd)
        
        return S
    
    @property
    def bedload_discharge(self):
        
        Qs = self._bedload_discharge * (1 - self.beta * self.distance__nd)
        
        return Qs
    
    @property
    def bedload_concentration(self):
        
        C = self._bedload_concentration * (1 - self.beta * self.distance__nd)
        
        return C
        
    @property
    def bankfull_discharge(self):
        
        return self._bankfull_discharge
        
        
    
    @property
    def bed_elevation(self):
        
        bed_elev = self.get_bed_elevation(self.time)
        
        return bed_elev
        
                
        
    @property
    def bankfull_channel_width(self):
        
        channel_width_bankfull = (
            (self.Cf / (self.coeff_Engelund * self.Shields_number**2.5 * sqrt(self.R))) *
            (self.bedload_discharge / (sqrt(self.g * self.grain_size) * self.grain_size)))
        
        return channel_width_bankfull
                
    @property         
    def bankfull_channel_depth(self):
        
        channel_depth_bankfull = (
            self.coeff_Engelund * (self.Shields_number**2  / sqrt(self.Cf)) *
            self._bankfull_discharge * self.grain_size / self.bedload_discharge)

        return channel_depth_bankfull
    