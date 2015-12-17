import numpy as np
from math import sqrt, log, exp
import FallVelocity

class DensityStratification(object):
    
    def __init__(self,
                 grain_size = 0.1,
                 flow_depth = 1.,
                 roughness_height = 0.002,
                 shear_velocity = 0.02,
                 skin_friction = 0.02,
                 reference_concentration = None,
                 rho_w = 1000.,
                 rho_s = 2650.,
                 g = 9.81,
                 kinematic_viscosity = 1e-6,
                 verbose = False,
                 save_output = True):
        
        
        self.verbose = verbose
        self.save_output = save_output

        self.grain_size = float(grain_size)
        self.flow_depth = float(flow_depth)
        self.kc = float(roughness_height)
        self.ustar = float(shear_velocity)
        self.nu = float(kinematic_viscosity)
        self.ustar_skin_friction = float(skin_friction)
        self.g = float(g)
        
        self.conc_ref = reference_concentration
        if reference_concentration:
            self.conc_ref = float(reference_concentration)


        self.kappa = 0.4 # von Karman's constant
        self.max_iterations = 200
        self.error_threshold = 0.001

        eta_ref = 0.05 # normalized reference bed elevation
        num_z_intervals = 50.
        self.dz = (1 - eta_ref) / num_z_intervals

        self.R = (float(rho_s) - float(rho_w)) / float(rho_w)

        self.eta = np.array(range(int(num_z_intervals)+1)) * self.dz + eta_ref
        self.Ri = np.zeros_like(self.eta)
        self.Fstrat = np.ones_like(self.eta)
        self.velocity_norm = np.zeros_like(self.eta)
        self.conc_norm = np.zeros_like(self.eta)
        self.int_c = np.zeros_like(self.eta)

        
        
    @property
    def gradient_Richardson_number(self):
        
        return self.Ri
        
    @property
    def normalized_velocity_no_stratification(self):
        
        return self.velocity_flat
    
    @property
    def normalized_concentration_no_stratification(self):
        
        return self.conc_norm
    
    @property
    def normalized_velocity_with_stratification(self):
        
        return self.velocity_norm
    
    @property
    def normalized_concentration_with_stratification(self):
        
        return self.conc_flat
        

    def compute_UC_normal(self):


        for i in range(1,len(self.eta)):

            ku1 = 1 / (self.kappa * self.eta[i-1] * self.Fstrat[i-1])
            ku2 = 1 / (self.kappa * self.eta[i] * self.Fstrat[i])
            kc1 = 1 / ((1.1 - self.eta[i-1]) * self.eta[i-1] * self.Fstrat[i-1])
            kc2 = 1 / ((1.1 - self.eta[i]) * self.eta[i] * self.Fstrat[i])

            self.velocity_norm[i] = self.velocity_norm[i-1] + (0.5 * (ku1 + ku2) * self.dz)
            self.int_c[i] = self.int_c[i-1] + (0.5 * (kc1 + kc2) * self.dz)

        self.conc_norm = np.exp(-self.int_c / (self.kappa * self.ustar_r))

        self.velocity_flat = 0.5 * self.dz * (self.velocity_norm[0] + self.velocity_norm[-1])
        self.conc_flat = 0.5 * self.dz * (self.conc_norm[0] + self.conc_norm[-1])
        self.qs_total = 0.5 * self.dz * (self.velocity_norm[0] * self.conc_norm[0] +
                                         self.velocity_norm[-1] * self.conc_norm[-1])

        self.velocity_flat = sum(self.dz * self.velocity_norm[1:]) + self.velocity_flat
        self.conc_flat = sum(self.dz * self.conc_norm[1:]) + self.conc_flat
        self.qs_total = sum(self.dz * self.velocity_norm[1:] * self.conc_norm[1:]) + self.qs_total

        self.Ri = (self.Ri_star * self.conc_norm * self.kappa * self.eta *
                   self.Fstrat / (self.ustar_r * (1.1 - self.eta)))
        X = 1.35 * self.Ri / (1 + (1.35 * self.Ri))
        self.Fstrat = 1 / (1 + (10. * X))

        self.Ri[-1] = 0
        self.Fstrat[-1] = 1



    def run(self):

        f = FallVelocity.FallVelocity(D = self.grain_size)
        f.run()

        Re = f.Reynolds_number
        Hr = self.flow_depth / self.kc
        unr = (1.0 / 0.4) * log(30.0 * 0.05 * Hr)
        vs = f.settling_velocity
        self.ustar_r = self.ustar / vs

        self.velocity_norm[0] = unr
        self.conc_norm[0] = 1
        self.Fstrat[0] = 1

        if self.conc_ref is None:

            Aa = 0.00000013

            # use the relation of Garcia and Parker (1991)
            Zgp = (self.ustar_skin_friction / vs) * Re**2./3
            self.conc_ref = Aa * Zgp**5. / (1 + (Aa / 0.3) * Zgp**5.)

        # Ri = gradient Richardson number
        self.Ri_star = self.R * 9.81 * self.flow_depth * self.conc_ref / self.ustar**2

        self.compute_UC_normal()

        # save initial conditions
        self.velocity_norm_init = self.velocity_norm.copy()
        self.conc_norm_init = self.conc_norm.copy()
        self.velocity_flat_init = self.velocity_flat
        self.conc_flat_init = self.conc_flat
        self.qs_total_init = self.qs_total


        flag_converged = False
        iteration = 0

        while not flag_converged and iteration <= self.max_iterations:

            iteration = iteration+1

            self.velocity_norm_o = self.velocity_norm.copy()
            self.conc_norm_o = self.conc_norm.copy()

            self.compute_UC_normal()

            error_vel = max(abs(2. * (self.velocity_norm - self.velocity_norm_o) /
                                     (self.velocity_norm + self.velocity_norm_o)))
            error_conc = max(abs(2. * (self.conc_norm - self.conc_norm_o) / 
                                      (self.conc_norm + self.conc_norm_o)))
            error = max((error_vel, error_conc))

            if error < self.error_threshold:
                flag_converged = True


        assert flag_converged, 'The calculation failed to converge.'
        
        
    def finalize():
        
        print 'finalize'
