import unittest
from AgDegNormalSub import AgDegNormalSub

class test_AgDegNormalSub(unittest.TestCase):
    
    def test_eta(self):
    
        f = AgDegNormalSub(friction = 'Manning',
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
							 verbose = False,
							 save_output = False)
                 
        f.run()
        self.assertEqual(f.bed_elevation[0], 102.26911259423767)
        self.assertEqual(f.bed_elevation[-2], 0.23724567562372906)
        
        
    def test_water_surface(self):
    
        f = AgDegNormalSub(friction = 'Manning',
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
							 verbose = False,
							 save_output = False)
        f.run()
        self.assertEqual(f.water_surface_elevation[0], 103.00768207498517)
        self.assertEqual(f.water_surface_elevation[-2], 1.100884370946358)
        
        
    def test_shear_stress(self):
    
        f = AgDegNormalSub(friction = 'Manning',
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
							 verbose = False,
							 save_output = False)
        f.run()
        self.assertEqual(f.bed_shear_stress[0], 0.21046460651748258)
        self.assertEqual(f.bed_shear_stress[-2], 0.14610015954984537)
        

if __name__ == '__main__':
    unittest.main()