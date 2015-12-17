import unittest
from Acronym import Acronym

class test_Acronym(unittest.TestCase):
    
    def test_percentfiner(self):
    
        default_GSD = [(40,100),
					   (35,99),
					   (30,83),
					   (25,72),
					   (20,41),
					   (15,22),
					   (10,5),
					   (5,0)]
					   
        f = Acronym(GSD = default_GSD,
                 ustar = 0.1,
                 rho_s = 2650.,
                 rho_w = 1000.,
                 verbose = False,
                 save_output = False)
        f.run()
        self.assertEqual(f.percent_finer['bedload'][2], 44.800790157927722)
        
    def test_characteristic_size(self):
    
        default_GSD = [(40,100),
					   (35,99),
					   (30,83),
					   (25,72),
					   (20,41),
					   (15,22),
					   (10,5),
					   (5,0)]
					   
        f = Acronym(GSD = default_GSD,
                 ustar = 0.1,
                 rho_s = 2650.,
                 rho_w = 1000.,
                 verbose = False,
                 save_output = False)
        f.run()
        self.assertEqual(f.characteristic_size['surface'][0], 16.93156268777064)
        self.assertEqual(f.characteristic_size['bedload'][3], 28.311241170777325)
        
    def test_stats(self):
    
        default_GSD = [(40,100),
					   (35,99),
					   (30,83),
					   (25,72),
					   (20,41),
					   (15,22),
					   (10,5),
					   (5,0)]
	
        f = Acronym(GSD = default_GSD,
                 ustar = 0.1,
                 rho_s = 2650.,
                 rho_w = 1000.,
                 verbose = False,
                 save_output = False)
        f.run()
        self.assertEqual(f.surface_statistics['sigma'], 0.5673872986902289)
        self.assertEqual(f.bedload_statistics['Dg'], 15.398111822119063)
        
        
if __name__ == '__main__':
    unittest.main()