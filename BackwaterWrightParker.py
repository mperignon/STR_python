import numpy as np
import json
from math import pow, exp, log, sqrt

class BackwaterWrightParker(object):
    """
    BackwaterWrightParker

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Dec 2015

    ----------------------------------------------

    Calculates backwater curves in sand-bed streams, including the effects of 
    both skin friction and form drag due to skin friction. Uses the hydraulic
    resistance formulation of Wright and Parker (2004) without the flow
    stratification correction

    """
        
    def __init__(self,
                 total_discharge = 1800.,
                 channel_width = 275.,
                 D50 = 0.35,
                 D90 = 1.,
                 bed_slope = 0.0001,
                 water_surface_elev__downstream = 20.,
                 reach_length = 200000.,
                 dx = 4000.,
                 rho_w = 1000.,
                 rho_s = 2650.,
                 g = 9.81,
                 save_output = True,
                 verbose = False):
        
        self.g = float(g)
        
        self.save_output = save_output
        self.verbose = verbose

        self.D50 = float(D50) / 1000.
        self.D90 = float(D90) / 1000.
        self.bed_slope = float(bed_slope)
        self.dx = dx
        
        self.discharge = float(total_discharge) / float(channel_width)
        self.R = (float(rho_s) - float(rho_w)) / float(rho_w)

        
        
        num_nodes = int(reach_length / self.dx + 1)
        self.num_iterations = 200
        self.K = 8.32

        # check that Froude number < 1
        self.Fr_number = self.discharge / (float(water_surface_elev__downstream) *
                                           sqrt(self.g * float(water_surface_elev__downstream)))
                                           
        assert self.Fr_number < 1, 'Froude number must be less than 1. Modify the water surface elevation.'

        # set initial bed
        self.x = np.arange(0, num_nodes*self.dx, self.dx)
        self.bed_elev = self.bed_slope * reach_length - self.bed_slope * self.x
        self.water_surface_elev = self.bed_elev.copy()
        self.water_surface_elev[-1] = float(water_surface_elev__downstream)
        self.water_depth = self.water_surface_elev - self.bed_elev
        
        self.water_slope = np.zeros_like(self.x)
        self.water_depth__sf = np.zeros_like(self.x)
        
        
    def run(self):

        # find normal depth
        self.water_depth__normal = self.Find_Normal_Depth()

        # calculate fluvial backwater
        self.Find_Fluvial_Backwater()

        
        
        
    def Find_Tau_Min(self, Fr):

        error_threshold = 0.001

        xtausmin = 0.4
        error = 9999.
        i = 0
        num_iterations = self.num_iterations

        while (error > error_threshold) and (i < num_iterations):

            i += 1

            Ft = xtausmin - 0.05 - 0.7 * pow(xtausmin, 0.8) * pow(Fr, (14./25))
            Ftp = 1 - 0.7 * 0.8 * pow(xtausmin, -0.2) * pow(Fr, (14./25))
            xtausmin_new = xtausmin - (Ft/Ftp)

            error = abs(2 * (xtausmin_new - xtausmin) / (xtausmin_new + xtausmin))

            xtausmin = xtausmin_new

        assert i <= num_iterations, 'Calculation of Tau min did not converge'

        return xtausmin


    def Find_Normal_Depth(self):

        error_threshold = 0.001
        error = 9999.
        i = 0
        num_iterations = self.num_iterations

        Hnorm = pow(pow(3. * self.D90, 1./3) * self.discharge**2 /
                    (self.K**2 * self.g * self.bed_slope), 0.3)

        while (error > error_threshold) and (i < num_iterations):

            i += 1

            Fr = self.discharge  / (sqrt(self.g) * pow(Hnorm, 1.5))
            tau = Hnorm * self.bed_slope / (self.R * self.D50)

            tau_min = self.Find_Tau_Min(Fr)

            H_threshold = self.R * self.D50 * tau_min / self.bed_slope
            fi = pow(self.discharge / (Hnorm * sqrt(self.g * Hnorm * self.bed_slope) *
                                       self.K * pow(Hnorm / (3. * self.D90), 1./6)), 1.5)

            FNt = Hnorm * fi
            FNtp = fi * (-1.5)

            if Hnorm >= H_threshold:

                FNt = (FNt - (self.R * self.D50 / self.bed_slope) *
                       (0.05 + 0.7 * pow(tau, 0.8) * pow(Fr, 14./25)))
                FNtp = (FNtp + self.R * self.D50 * 0.7 * pow(tau, 0.8) *
                        pow(Fr, 14./25) / (25. * Hnorm * self.bed_slope))

            else:

                FNt = FNt - Hnorm
                FNtp = FNtp - 1


            Hnorm_new = Hnorm - FNt / FNtp

            error = abs(2 * (Hnorm_new - Hnorm) / (Hnorm_new + Hnorm))

            Hnorm = Hnorm_new

        assert i <= num_iterations, 'Calculation of Normal Depth didnot converge'

        return Hnorm
    
    @property
    def depth(self):
        return self.water_depth
    
    @property
    def water_surface_elevation(self):
        return self.water_surface_elev
    
    @property
    def bed_elevation(self):
        return self.bed_elev
    
    @property
    def downstream_water_surface_elevation(self):
        return self.water_surface_elev__ds
    
    @property
    def Froude_number(self):
        return self.Fr_number



    def Compute_Sf_and_Hs(self, H):

        Fr_loc = self.discharge / (sqrt(self.g) * pow(H, 3./2))

        Snom = pow(Fr_loc / (self.K * pow(H / (3. * self.D90), 1./6)), 2)
        tau_nom = H * Snom / (self.R * self.D50)

        fi_loc = 0.9

        error_threshold = 0.001
        error = 9999.
        i = 0
        num_iterations = 20

        while (error > error_threshold) and (i < num_iterations):

            i += 1

            tau_min = self.Find_Tau_Min(Fr_loc)

            ratio_c = pow(tau_nom / tau_min, 0.75)

            a1 = 0.7 * pow(tau_nom, 0.8) * pow(Fr_loc, 14./25)

            if fi_loc <= ratio_c:

                a2 = tau_nom / pow(fi_loc, 1./3) - 0.05
                Ff = fi_loc - pow(a2/a1, -15./16)
                Ffp = 1 - (15./16) * pow(a2/a1, -31./16) * pow(fi_loc, -4./3) * tau_nom / (3. * a1)

            else:

                Ff = fi_loc - 1
                Ffp = 1.

            fi_loc_new = fi_loc - (Ff/Ffp)

            if fi_loc_new < 0:
                fi_loc_new = 1.02

            error = abs(2. * (fi_loc_new - fi_loc) / (fi_loc_new + fi_loc))

            fi_loc = fi_loc_new


        assert i <= num_iterations, 'Calculation of Sf and Hs did not converge'

        Hs = fi_loc * H
        Sf = pow(fi_loc, -4./3) * Snom

        return Hs, Sf
                        

    def Find_Fluvial_Backwater(self):

        self.water_depth__sf[-1], self.water_slope[-1] = self.Compute_Sf_and_Hs(self.water_depth[-1])
        
        for i in range(len(self.x)-1):

            fr2dn = self.discharge**2 / (self.g * pow(self.water_depth[-1-i],3))
            fnp = (self.bed_elev[-2-i] - self.bed_elev[-1-i] -
                   self.water_slope[-1-i]*self.dx) / (1 - fr2dn)

            H_pred = self.water_depth[-1-i] - fnp
            Hs_pred, S_pred = self.Compute_Sf_and_Hs(H_pred)

            fr2up = self.discharge**2 / (self.g * H_pred**3)
            fn = (self.bed_elev[-2-i] - self.bed_elev[-1-i] -
                  S_pred*self.dx) / (1 - fr2up)

            self.water_depth[-2-i] = self.water_depth[-1-i] - 0.5 * (fnp + fn)

            self.water_depth__sf[-2-i], self.water_slope[-2-i] = self.Compute_Sf_and_Hs(self.water_depth[-2-i])


        self.water_surface_elev = self.bed_elev + self.water_depth

        
    ##########################################
    ############## Finalize ##################
    ##########################################   
    def finalize(self):
        
        if self.save_output:
            
            output_dict = {
                'Water_depth_normal' : self.water_depth__normal,
                'Froude_number' : self.Fr_number
            }
            
            with open('output/BackwaterWP_normal.json', 'w') as f:
                json.dump(output_dict, f, indent=4)
                


            fields = ['Position', 'Bed_elevation', 'Water_surface_elevation',
                      'Depth', 'Depth_from_skin_friction']

            header = ', '.join(fields)

            dtype = [(i, float) for i in fields]
            self.data = np.empty(len(self.x), dtype = dtype)

            self.data['Position'] = self.x
            self.data['Bed_elevation'] = self.bed_elev
            self.data['Water_surface_elevation'] = self.water_surface_elev
            self.data['Depth'] = self.water_depth
            self.data['Depth_from_skin_friction'] = self.water_depth__sf

            np.savetxt('output/BackwaterWG_profile.csv', self.data,
                       header = header,
                       delimiter = ",",
                       fmt = '%10.5f',
                       comments = '')

            if self.verbose:
                print 'Saved data to file'