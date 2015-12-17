import unittest
from AgDegNormal import AgDegNormal

class test_AgDegNormal(unittest.TestCase):
    
    def test_eta(self):
    
        f = AgDegNormal(friction = 'Manning',
						 Qf = 70,
						 If = 0.03,
						 B = 25,
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
						 verbose = False,
						 save_output = False)
                 
        f.run()
        self.assertEqual(f.bed_elevation[0], 141.47327348563991)
        self.assertEqual(f.bed_elevation[-2], 1.4114131765610818)
        
        
    def test_water_surface(self):
    
        f = AgDegNormal(friction = 'Manning',
						 Qf = 70,
						 If = 0.03,
						 B = 25,
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
						 verbose = False,
						 save_output = False)
        f.run()
        self.assertEqual(f.water_surface_elevation[0], 142.21027987395911)
        self.assertEqual(f.water_surface_elevation[-2], 2.1498495728966525)
        
        
    def test_shear_stress(self):
    
        f = AgDegNormal(friction = 'Manning',
						 Qf = 70,
						 If = 0.03,
						 B = 25,
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
						 verbose = False,
						 save_output = False)
        f.run()
        self.assertEqual(f.bed_shear_stress[0], 0.21150760344400008)
        self.assertEqual(f.bed_shear_stress[-2], 0.21055312258898984)
        

if __name__ == '__main__':
    unittest.main()