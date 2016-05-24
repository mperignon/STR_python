import unittest
from SteadyStateAg import SteadyStateAg

class test_SteadyStateAg(unittest.TestCase):
    
    def test_bedload_discharge(self):
        f = SteadyStateAg(time = 10000.,
                         channel_length = 400000.,
                         floodplain_width = 8000.,
                         channel_sinuosity = 1.5,
                         bed_porosity = 0.4,
                         flood_intermittence = 0.1,
                         bankfull_discharge = 5000.,
                         annual_bedload_discharge = 12.,
                         grain_size = 0.25,
                         fraction_washload_deposition = 2.5,
                         sea_level_elevation = 0.,
                         sea_level_change = 0.010,
                         Shields_number = 1.86,
                         coeff_Engelund = 0.05,
                         coeff_Chezy = 25.,
                         density_water = 1000.,
                         density_sediment = 2650.,
                         gravity = 9.81,
                         verbose = False,
                         save_output = False)
        f.run()
                                    
        self.assertEqual(f.bedload_discharge[0], 1.4349322783711222)
        self.assertEqual(f.bedload_discharge[-2], 0.33399757031800209)
        
        
    def test_bed_elev_change(self):
        f = SteadyStateAg(time = 10000.,
                         channel_length = 400000.,
                         floodplain_width = 8000.,
                         channel_sinuosity = 1.5,
                         bed_porosity = 0.4,
                         flood_intermittence = 0.1,
                         bankfull_discharge = 5000.,
                         annual_bedload_discharge = 12.,
                         grain_size = 0.25,
                         fraction_washload_deposition = 2.5,
                         sea_level_elevation = 0.,
                         sea_level_change = 0.010,
                         Shields_number = 1.86,
                         coeff_Engelund = 0.05,
                         coeff_Chezy = 25.,
                         density_water = 1000.,
                         density_sediment = 2650.,
                         gravity = 9.81,
                         verbose = False,
                         save_output = False)
        f.run()
                                    
        self.assertEqual(f.bed_elevation_change[0], 48.56992279625706)
        self.assertEqual(f.bed_elevation_change[-2], 0.86587913483422263)
        
        
    def test_bed_elev(self):
        f = SteadyStateAg(time = 10000.,
                         channel_length = 400000.,
                         floodplain_width = 8000.,
                         channel_sinuosity = 1.5,
                         bed_porosity = 0.4,
                         flood_intermittence = 0.1,
                         bankfull_discharge = 5000.,
                         annual_bedload_discharge = 12.,
                         grain_size = 0.25,
                         fraction_washload_deposition = 2.5,
                         sea_level_elevation = 0.,
                         sea_level_change = 0.010,
                         Shields_number = 1.86,
                         coeff_Engelund = 0.05,
                         coeff_Chezy = 25.,
                         density_water = 1000.,
                         density_sediment = 2650.,
                         gravity = 9.81,
                         verbose = False,
                         save_output = False)
        f.run()
                                    
        self.assertEqual(f.bed_elevation[0][0], 148.56992279625706)
        self.assertEqual(f.bed_elevation[0][-2], 100.86587913483422)
        
        
    def test_bankfull_width(self):
        f = SteadyStateAg(time = 10000.,
                         channel_length = 400000.,
                         floodplain_width = 8000.,
                         channel_sinuosity = 1.5,
                         bed_porosity = 0.4,
                         flood_intermittence = 0.1,
                         bankfull_discharge = 5000.,
                         annual_bedload_discharge = 12.,
                         grain_size = 0.25,
                         fraction_washload_deposition = 2.5,
                         sea_level_elevation = 0.,
                         sea_level_change = 0.010,
                         Shields_number = 1.86,
                         coeff_Engelund = 0.05,
                         coeff_Chezy = 25.,
                         density_water = 1000.,
                         density_sediment = 2650.,
                         gravity = 9.81,
                         verbose = False,
                         save_output = False)
        f.run()
                                    
        self.assertEqual(f.bankfull_channel_width[0], 611.94520196967449)
        self.assertEqual(f.bankfull_channel_width[-2], 142.43753082036991)
        
        
    def test_bankfull_depth(self):
        f = SteadyStateAg(time = 10000.,
                         channel_length = 400000.,
                         floodplain_width = 8000.,
                         channel_sinuosity = 1.5,
                         bed_porosity = 0.4,
                         flood_intermittence = 0.1,
                         bankfull_discharge = 5000.,
                         annual_bedload_discharge = 12.,
                         grain_size = 0.25,
                         fraction_washload_deposition = 2.5,
                         sea_level_elevation = 0.,
                         sea_level_change = 0.010,
                         Shields_number = 1.86,
                         coeff_Engelund = 0.05,
                         coeff_Chezy = 25.,
                         density_water = 1000.,
                         density_sediment = 2650.,
                         gravity = 9.81,
                         verbose = False,
                         save_output = False)
        f.run()
                                    
        self.assertEqual(f.bankfull_channel_depth[0], 3.7671638456250012)
        self.assertEqual(f.bankfull_channel_depth[-2], 16.184623723020671)
 
        
if __name__ == '__main__':
    unittest.main()