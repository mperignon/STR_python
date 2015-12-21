import unittest
from BackwaterWrightParker import BackwaterWrightParker

class test_BackwaterWrightParker(unittest.TestCase):
    
    def test_depth(self):
        f = BackwaterWrightParker(total_discharge = 1800.,
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
                                 save_output = False,
                                 verbose = False)
        f.run()
        self.assertEqual(f.water_depth[0], 6.7941284285573964)
        
    def test_depth_sf(self):
        f = BackwaterWrightParker(total_discharge = 1800.,
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
                                 save_output = False,
                                 verbose = False)
        f.run()
        self.assertEqual(f.water_depth__sf[0], 1.6934703386292611)
        
    def test_Fr(self):
        f = BackwaterWrightParker(total_discharge = 1800.,
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
                                 save_output = False,
                                 verbose = False)
        f.run()
        self.assertEqual(f.Froude_number, 0.02336470564789742)

        
if __name__ == '__main__':
    unittest.main()