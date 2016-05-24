import unittest
from AgDegNormalFault import AgDegNormalFault

class test_AgDegNormalFault(unittest.TestCase):
    
    def test_eta(self):
    
        f = AgDegNormalFault(friction = 'Manning',
							 Qf = 70,
							 If = 0.03,
							 L = 10000,
							 dx = 100,
							 S = 0.008,
							 Gtf = 700000,
							 D = 30,
							 porosity = 0.35,
							 kc = 75,
							 phi = 1,
							 alpha_r = 8.1,
							 alpha_t = 8.,
							 tau_c = 0.047,
							 nt = 1.5,
							 Cz = 22,
							 alpha_u = 0.5,
							 dt = 0.01,
							 totalT = 120,
							 print_dt = 120,
							 g = 9.81,
							 rho_w = 1000,
							 rho_s = 2650,
							 subsidence = 5,
							 channel_width = 25,
							 depositional_width = 25*60,
							 washload_fraction = 1.5,
							 sinuosity = 1.5,
							 fault_position = 0.5,
							 fault_height = 15,
							 fault_time = 1,
							 verbose = False,
							 save_output = False)
                 
        f.run()
        self.assertEqual(f.bed_elevation[0], 126.64142774887191)
        self.assertEqual(f.bed_elevation[-2], -13.585970469445467)
        
        
    def test_water_surface(self):
    
        f = AgDegNormalFault(friction = 'Manning',
							 Qf = 70,
							 If = 0.03,
							 L = 10000,
							 dx = 100,
							 S = 0.008,
							 Gtf = 700000,
							 D = 30,
							 porosity = 0.35,
							 kc = 75,
							 phi = 1,
							 alpha_r = 8.1,
							 alpha_t = 8.,
							 tau_c = 0.047,
							 nt = 1.5,
							 Cz = 22,
							 alpha_u = 0.5,
							 dt = 0.01,
							 totalT = 120,
							 print_dt = 120,
							 g = 9.81,
							 rho_w = 1000,
							 rho_s = 2650,
							 subsidence = 5,
							 channel_width = 25,
							 depositional_width = 25*60,
							 washload_fraction = 1.5,
							 sinuosity = 1.5,
							 fault_position = 0.5,
							 fault_height = 15,
							 fault_time = 1,
							 verbose = False,
							 save_output = False)
        f.run()
        self.assertEqual(f.water_surface_elevation[0], 127.37842780835773)
        self.assertEqual(f.water_surface_elevation[-2], -12.847944287798242)
        
        
    def test_shear_stress(self):
    
        f = AgDegNormalFault(friction = 'Manning',
							 Qf = 70,
							 If = 0.03,
							 L = 10000,
							 dx = 100,
							 S = 0.008,
							 Gtf = 700000,
							 D = 30,
							 porosity = 0.35,
							 kc = 75,
							 phi = 1,
							 alpha_r = 8.1,
							 alpha_t = 8.,
							 tau_c = 0.047,
							 nt = 1.5,
							 Cz = 22,
							 alpha_u = 0.5,
							 dt = 0.01,
							 totalT = 120,
							 print_dt = 120,
							 g = 9.81,
							 rho_w = 1000,
							 rho_s = 2650,
							 subsidence = 5,
							 channel_width = 25,
							 depositional_width = 25*60,
							 washload_fraction = 1.5,
							 sinuosity = 1.5,
							 fault_position = 0.5,
							 fault_height = 15,
							 fault_time = 1,
							 verbose = False,
							 save_output = False)
        f.run()
        self.assertEqual(f.bed_shear_stress[0], 0.21151184144832058)
        self.assertEqual(f.bed_shear_stress[-2], 0.21082629624709234)
        

if __name__ == '__main__':
    unittest.main()