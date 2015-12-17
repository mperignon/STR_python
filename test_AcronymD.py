import unittest
from Acronym_D import Acronym_D

class test_Acronym_D(unittest.TestCase):
    
    def test_percentfiner(self):
    
		default_GSD = [(32,100),(16,50),(8,0)]
					   
					   
		default_flow = [(80,0),
						(40, 20),
						(30, 60),
						(10, 100)]

		f = Acronym_D(GSD = default_GSD,
					 flow_curve = default_flow,
					 channel_slope = 0.005,
					 channel_width = 12.,
					 mannings_n = 2.,
					 rho_s = 2650.,
					 rho_w = 1000.,
					 alpha_r = 8.1,
					 verbose = False,
					 save_output = False)
					 
		f.run()
		self.assertEqual(f.percent_finer['bedload'][1], 52.424002448886377)
		
		
    def test_characteristic_size(self):
    
		default_GSD = [(32,100),(16,50),(8,0)]
					   
					   
		default_flow = [(80,0),
						(40, 20),
						(30, 60),
						(10, 100)]

		f = Acronym_D(GSD = default_GSD,
					 flow_curve = default_flow,
					 channel_slope = 0.005,
					 channel_width = 12.,
					 mannings_n = 2.,
					 rho_s = 2650.,
					 rho_w = 1000.,
					 alpha_r = 8.1,
					 verbose = False,
					 save_output = False)
					 
		f.run()
		self.assertEqual(f.characteristic_size['surface'][0], 12.125732532083186)
		self.assertEqual(f.characteristic_size['bedload'][3], 27.661547759243085)
        
    def test_stats(self):
    
		default_GSD = [(32,100),(16,50),(8,0)]
					   
					   
		default_flow = [(80,0),
						(40, 20),
						(30, 60),
						(10, 100)]

		f = Acronym_D(GSD = default_GSD,
					 flow_curve = default_flow,
					 channel_slope = 0.005,
					 channel_width = 12.,
					 mannings_n = 2.,
					 rho_s = 2650.,
					 rho_w = 1000.,
					 alpha_r = 8.1,
					 verbose = False,
					 save_output = False)
		f.run()
		self.assertEqual(f.surface_statistics['sigma_g'], 1.4142135623730951)
		self.assertEqual(f.bedload_statistics['Dg'], 15.733415362194647)
        
    def test_averages(self):
    
		default_GSD = [(32,100),(16,50),(8,0)]
					   
					   
		default_flow = [(80,0),
						(40, 20),
						(30, 60),
						(10, 100)]

		f = Acronym_D(GSD = default_GSD,
					 flow_curve = default_flow,
					 channel_slope = 0.005,
					 channel_width = 12.,
					 mannings_n = 2.,
					 rho_s = 2650.,
					 rho_w = 1000.,
					 alpha_r = 8.1,
					 verbose = False,
					 save_output = False)
		f.run()
		self.assertEqual(f.averages['H'][0], 0.99594590233287528)
		self.assertEqual(f.averages['ustar'][0], 0.20904528549189427)


        
if __name__ == '__main__':
    unittest.main()