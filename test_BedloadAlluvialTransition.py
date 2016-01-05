import unittest
from BedrockAlluvialTransition import BedrockAlluvialTransition

class test_BedrockAlluvialTransition(unittest.TestCase):
    
    def test_depth(self):
        f = BedrockAlluvialTransition(flood_discharge = 6.,
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
                                     total_time = 100,
                                     print_dt = 100,
                                     critical_shear_stress = 0.,
                                     rho_w = 1000.,
                                     rho_s = 2650.,
                                     verbose = False,
                                     save_output = False)
        f.run()
        self.assertEqual(f.water_depth[1], 4.3002211211040553)


    def test_eta(self):
        f = BedrockAlluvialTransition(flood_discharge = 6.,
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
                                     total_time = 100,
                                     print_dt = 100,
                                     critical_shear_stress = 0.,
                                     rho_w = 1000.,
                                     rho_s = 2650.,
                                     verbose = False,
                                     save_output = False)
        f.run()
        self.assertEqual(f.bed_elevation[1], 13.579561099765957)
   
        
    def test_x(self):
        f = BedrockAlluvialTransition(flood_discharge = 6.,
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
                                     total_time = 100,
                                     print_dt = 100,
                                     critical_shear_stress = 0.,
                                     rho_w = 1000.,
                                     rho_s = 2650.,
                                     verbose = False,
                                     save_output = False)
        f.run()
        self.assertEqual(f.distance[1], -24322.229449656246)        
        
if __name__ == '__main__':
    unittest.main()