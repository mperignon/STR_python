import unittest
from AgDegBW import AgDegBW

class test_AgDegBW(unittest.TestCase):
    
    def test_eta(self):
    
        f = AgDegBW(friction = "Manning",
                     Qf = 400,
                     If = 0.1,
                     S = 0.0002,
                     initdepth = 15,
                     B = 60,
                     L = 50000,
                     dx = 833,
                     Gtf = 187059,
                     D90 = 0.6,
                     porosity = 0.4,
                     kc = 75,
                     phi = 1,
                     alpha_r = 8.1,
                     alpha_t = 8.,
                     tau_c = 0.047,
                     n = 1.5,
                     Cz = 22,
                     alpha_u = 1,
                     dt = 0.1,
                     totalT = 120,
                     print_dt = 120,
                     g = 9.81,
                     rho_w = 1000,
                     rho_s = 2650,
                     verbose = False,
                     save_output = False)
                 
        f.run()
        self.assertEqual(f.bed_elevation[0], 16.150678907133734)
        self.assertEqual(f.bed_elevation[-2], 0.23969483388856572)
        
        
    def test_water_surface(self):
    
        f = AgDegBW(friction = "Manning",
                     Qf = 400,
                     If = 0.1,
                     S = 0.0002,
                     initdepth = 15,
                     B = 60,
                     L = 50000,
                     dx = 833,
                     Gtf = 187059,
                     D90 = 0.6,
                     porosity = 0.4,
                     kc = 75,
                     phi = 1,
                     alpha_r = 8.1,
                     alpha_t = 8.,
                     tau_c = 0.047,
                     n = 1.5,
                     Cz = 22,
                     alpha_u = 1,
                     dt = 0.1,
                     totalT = 120,
                     print_dt = 120,
                     g = 9.81,
                     rho_w = 1000,
                     rho_s = 2650,
                     verbose = False,
                     save_output = False)
        f.run()
        self.assertEqual(f.water_surface_elevation[0], 20.619136144306204)
        self.assertEqual(f.water_surface_elevation[-2], 15.002769439775431)
        
        
    def test_shear_stress(self):
    
        f = AgDegBW(friction = "Manning",
                     Qf = 400,
                     If = 0.1,
                     S = 0.0002,
                     initdepth = 15,
                     B = 60,
                     L = 50000,
                     dx = 833,
                     Gtf = 187059,
                     D90 = 0.6,
                     porosity = 0.4,
                     kc = 75,
                     phi = 1,
                     alpha_r = 8.1,
                     alpha_t = 8.,
                     tau_c = 0.047,
                     n = 1.5,
                     Cz = 22,
                     alpha_u = 1,
                     dt = 0.1,
                     totalT = 120,
                     print_dt = 120,
                     g = 9.81,
                     rho_w = 1000,
                     rho_s = 2650,
                     verbose = False,
                     save_output = False)
        f.run()
        self.assertEqual(f.bed_shear_stress[0], 0.89439215587050513)
        self.assertEqual(f.bed_shear_stress[-2], 0.055014859454631616)
        

if __name__ == '__main__':
    unittest.main()