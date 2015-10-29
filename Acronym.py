import numpy as np
import matplotlib.pyplot as plt

import GSD_calculator as gsd

from math import sqrt, log, exp
from shared_tools import po, oo, so
import shared_tools as tools

default_GSD = [(40,100),
               (35,99),
               (30,83),
               (25,72),
               (20,41),
               (15,22),
               (10,5),
               (5,0)]



class Acronym(object):


    def __init__(self,
                 GSD = default_GSD,
                 ustar = 0.1,
                 rho_s = 2650.,
                 rho_w = 1000.):

        
        self.ustar = float(ustar)
        self.GSD = GSD

        self.R = (float(rho_s) - float(rho_w))/ float(rho_w)
        self.g = 9.81

        self._percent_finer = None
        self._characteristic_size = None
        self.tau_star_surface = None
        self._surface_statistics = None
        self._bedload_statistics = None




    def run(self):

        self.GSD_edges, self.GSD_surface = self.process_surface_GSD(GSD)

        self.pf_surface = self.GSD_edges['fractions']*100
        
        
        if len(self.GSD_surface) > 1:
            
            surface_vals = self.grain_statistics(self.GSD_edges,
                                                 self.GSD_surface,
                                                 self.pf_surface)
            
            self.c_sizes_surface = surface_vals[0]
            self.D_g_surface = surface_vals[1]
            self.sigma_g_surface = surface_vals[2]
            self.sigma_surface = surface_vals[3]
            
            

        calc_bedload_vals = self.calculate_bedload(self.GSD_surface,
                                                   self.D_g_surface,
                                                   self.sigma_surface)
        
        self.tau_star_surface = calc_bedload_vals[0]
        self.GSD_bedload = calc_bedload_vals[1]
        self.q_bedload = calc_bedload_vals[2]
        self.q_total_bedload = calc_bedload_vals[3]


        self.pf_bedload = self.get_percent_finer(self.GSD_bedload)
        
        if len(self.GSD_bedload) > 1:
            
            bedload_vals = self.grain_statistics(self.GSD_edges,
                                                 self.GSD_bedload,
                                                 self.pf_bedload)
            
            self.c_sizes_bedload = bedload_vals[0]
            self.D_g_bedload = bedload_vals[1]
            self.sigma_g_bedload = bedload_vals[2]
            self.sigma_bedload = bedload_vals[3]
        

        
    @property
    def grain_size_distribution(self):
        
        return self.GSD
    
    @grain_size_distribution.setter
    def grain_size_distribution(self, new_GSD):
        
        self.GSD = new_GSD
        
        

    @property
    def percent_finer(self):
        
        titles = ('diameter', 'surface', 'bedload')
        data = (self.GSD_edges['diameters'], self.pf_surface, self.pf_bedload)
        self._percent_finer = self.create_output_array(titles, data)
        
        return self._percent_finer
    
    
    @property
    def characteristic_size(self):
        
        if len(self.GSD_edges) > 2:
            titles = ('percentile', 'surface', 'bedload')
            data = (np.array((30, 50, 70, 90)), self.c_sizes_surface, self.c_sizes_bedload)
            self._characteristic_size = self.create_output_array(titles, data)
        
        return self._characteristic_size
    
    
    @property
    def surface_statistics(self):
        
        self._surface_statistics = {'Dg': self.D_g_surface,
                                   'sigma': self.sigma_surface,
                                   'sigma_g': self.sigma_g_surface}
        
        return self._surface_statistics
    
    
    @property
    def bedload_statistics(self):
        
        self._bedload_statistics = {'Dg': self.D_g_bedload,
                                   'sigma': self.sigma_bedload,
                                   'sigma_g': self.sigma_g_bedload}
        
        return self._bedload_statistics
    
    
    @property
    def Shields_stress(self):
        
        return self.tau_star_surface
            

    def get_bedload_discharge(self, GSD, Omega, phi_star, D_g):

        if len(GSD) == 1:

            phi = Omega * phi_star
            q_total = self.GG(phi) * 0.00218 * self.ustar**3 / (self.R * self.g)
            q = q_total

        else:

            q = np.zeros((len(GSD),))

            for i in range(len(q)):

                zs = GSD['diameters'][i] / D_g
                phi = Omega * phi_star * zs**(-0.0951)
                q[i] = GSD['fractions'][i] * self.GG(phi) * \
                                0.00218 * self.ustar**3 / (self.R * self.g)

            q_total = sum(q)

        return q, q_total


    def get_Omega(self, phi_star, sigma):

        if phi_star > po[-2]:    
            omega_o = oo[-2] + ((oo[-1] - oo[-2])/(po[-1] - po[-2])) * (phi_star - po[-2])
            sigma_o = so[-2] + ((so[-1] - so[-2])/(po[-1] - po[-2])) * (phi_star - po[-2])

        elif phi_star <= 0.7639:    
            omega_o = 1.011
            sigma_o = 0.8157

        else:

            ind_po = np.where(phi_star > po)[0][-1]

            xm1 = po[ind_po-1]
            x = po[ind_po];
            xp1 = po[ind_po+1];
            xp2 = po[ind_po+2];

            ym1 = oo[ind_po-1];
            y = oo[ind_po];
            yp1 = oo[ind_po+1];
            yp2 = oo[ind_po+2];

            omega_o = tools.CubicInterpolant(ym1, y, yp1, yp2,
                                             xm1, x, xp1, xp2,
                                             phi_star);

            ym1 = so[ind_po-1];
            y = so[ind_po];
            yp1 = so[ind_po+1];
            yp2 = so[ind_po+2];

            sigma_o = tools.CubicInterpolant(ym1, y, yp1, yp2,
                                             xm1, x, xp1, xp2,
                                             phi_star);


        Omega = 1.0 + (sigma / sigma_o) * (omega_o - 1.0)

        return Omega


    def grain_statistics(self, GSD, centerpts, pf):

        psi_bar = sum(centerpts['psi'] * centerpts['fractions'])
        D_g = 2 ** psi_bar

        sigma = sqrt(sum((centerpts['psi'] - psi_bar)**2 * \
            centerpts['fractions']))
        sigma_g = 2 ** sigma

        i = 0
        psi = np.zeros((len(GSD),))
        
        for x in range(30,100,20):

            while x > GSD['fractions'][i]*100:
                i += 1

            psi[i-1] = log(GSD['diameters'][i-1], 2)
            psi[i] = log(GSD['diameters'][i], 2)
            psix = psi[i-1] + ((psi[i] - psi[i-1])/(pf[i] - pf[i-1])) * \
                    (x - pf[i-1])
            if x == 30:
                Dx30 = 2**psix
            if x == 50:
                Dx50 = 2**psix
            if x == 70:
                Dx70 = 2**psix
            else:
                Dx90 = 2**psix 

        return ((Dx30, Dx50, Dx70, Dx90), D_g, sigma_g, sigma)



    def GG(self, phi):

        Mo = 14.2

        if phi < 1:
            G = phi ** Mo
        elif phi <= 1.59 and phi >= 1:
            G = exp(Mo * (phi - 1.0) - 9.28 * (phi - 1.0)**2)
        else:
            G = 5474. * (1.0 - (0.853/phi))**4.5

        return G



    def create_center_array_from_GSD_edges(self, datapts):

        psi_b = [log(i,2) for i in datapts['diameters']]
        

        dtype = [('psi', float), ('diameters', float), ('fractions', float)]
        centerpts = np.zeros((len(datapts)-1), dtype=dtype)

        for i in range(len(centerpts)):    
            centerpts['psi'][i] = (psi_b[i+1] + psi_b[i]) / 2
            centerpts['diameters'][i] = sqrt(datapts['diameters'][i] * 
                                             datapts['diameters'][i+1])
            centerpts['fractions'][i] = (datapts['fractions'][i+1] - 
                                         datapts['fractions'][i])

                
        return centerpts



    def get_percent_finer(self, GSD):

        # percent finer
        pf = np.zeros((len(GSD)+1,))
        pf[-1] = 100

        for i in range(-2,-len(GSD)-1,-1):
            pf[i] = pf[i+1] - 100 * GSD['fractions'][i+1]

        return pf


    def calculate_bedload(self, GSD, D_g, sigma):
        
        taursgo = 0.0386
        lowphi = 0.1
        hiphi = 2320.

        tau_star = self.ustar**2 / (self.R * self.g * D_g/1000)
        phi_star = tau_star / taursgo

        assert phi_star > lowphi, "Friction velocity is too low."
        assert phi_star < hiphi, "Friction velocity is too high."

        Omega = self.get_Omega(phi_star, sigma)

        q, q_total = self.get_bedload_discharge(GSD,
                                                Omega,
                                                phi_star,
                                                D_g)
        

        titles = ('psi', 'diameters', 'fractions')
        
        data = (GSD['psi'], GSD['diameters'], q / q_total)
        GSD_bedload = self.create_output_array(titles, data)

        return (tau_star, GSD_bedload, q, q_total)


    def process_surface_GSD(self, GSD):

        assert len(GSD) >= 2, 'Too few grain sizes. At least two size bounds required.'

        GSD_edges = gsd.pre_process_data(GSD)
        GSD_surface = self.create_center_array_from_GSD_edges(GSD_edges)

        # check for sand. If have D<2 with fraction > 0.0, break
        if min(GSD_edges['diameters']) <= 2:

            error_msg = "This program cannot be used with sand-sized sediment"
            max_sand = max(GSD_edges[GSD_edges['diameters'] <= 2]['fractions'])
            
            assert max_sand <= 0.0, error_msg

        return GSD_edges, GSD_surface


    def create_output_array(self, titles, data):

        dtype = [(titles[0], float), (titles[1], float), (titles[2], float)]
        out_array = np.zeros((len(data[0])), dtype=dtype)

        out_array[titles[0]] = data[0]
        out_array[titles[1]] = data[1]
        out_array[titles[2]] = data[2]

        return out_array

