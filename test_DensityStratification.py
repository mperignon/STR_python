import unittest
from DensityStratification import DensityStratification

class test_DensityStratification(unittest.TestCase):
    
    def test_conc(self):
        f = DensityStratification(grain_size = 0.1,
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
                                     save_output = False)
        f.run()
        self.assertEqual(f.conc_norm_init[-1], 0.0009073334043798784)
        
    def test_velocity(self):
        f = DensityStratification(grain_size = 0.1,
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
                                     save_output = False)
        f.run()
        self.assertEqual(f.velocity_norm[0], 16.550183016325889)

        
if __name__ == '__main__':
    unittest.main()