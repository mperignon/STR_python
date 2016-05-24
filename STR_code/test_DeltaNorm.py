import unittest
from DeltaBW import DeltaBW

class test_DeltaBW(unittest.TestCase):
    
    def test_mass_balance(self):
        f = DeltaBW(friction = "Chezy",
                        flood_discharge = 6.,
                        flood_intermittency = 1.,
                        bedload_discharge_input = 0.001,
                        grain_size = 0.5,
                        coeff_Chezy = 15.,
                        exponent_load_relation = 2.5,
                        critical_Shields_stress = 0.,
                        foreset_elev__top__init = 3.,
                        foreset_elev__base__init = 0.,
                        bed_slope__init = 0.00025,
                        basin_slope = 0.,
                        domain_length = 10000.,
                        domain_length_max = 500000.,
                        foreset_slope = 0.2,
                        bed_porosity = 0.4,
                        coeff_Manning_ks = 0.3,
                        coeff_bedload = 7.2,
                        coeff_Manning = 8.1,
                        dx = 500.,
                        dt__days = 0.182625,
                        num_iterations = 10000,
                        print_step = 5000,
                        density_sediment = 2650.,
                        density_water = 1000.,
                        basin_water_surface_elevation = 8.5,
                        verbose = False,
                        save_output = False)
        f.run()
        f.finalize()
        self.assertEqual(f.sed_volume__init, 26338.5)
        self.assertEqual(f.sed_volume__final, 190680.45389578832)
        self.assertEqual(f.sed_volume__feed, 157788.00000000003)
        self.assertEqual(f.mass_balance_error, 0.034371398650908352)
        
    def test_stats(self):
        f = DeltaBW(friction = "Chezy",
                        flood_discharge = 6.,
                        flood_intermittency = 1.,
                        bedload_discharge_input = 0.001,
                        grain_size = 0.5,
                        coeff_Chezy = 15.,
                        exponent_load_relation = 2.5,
                        critical_Shields_stress = 0.,
                        foreset_elev__top__init = 3.,
                        foreset_elev__base__init = 0.,
                        bed_slope__init = 0.00025,
                        basin_slope = 0.,
                        domain_length = 10000.,
                        domain_length_max = 500000.,
                        foreset_slope = 0.2,
                        bed_porosity = 0.4,
                        coeff_Manning_ks = 0.3,
                        coeff_bedload = 7.2,
                        coeff_Manning = 8.1,
                        dx = 500.,
                        dt__days = 0.182625,
                        num_iterations = 10000,
                        print_step = 5000,
                        density_sediment = 2650.,
                        density_water = 1000.,
                        basin_water_surface_elevation = 8.5,
                        verbose = False,
                        save_output = False)
        f.run()
        self.assertEqual(f.domain_length, 34481.071182531079)
        self.assertEqual(f.bed_length, 34500.977813788544)
        self.assertEqual(f.eta_upstream, 15.154555352651583)
        self.assertEqual(f.eta[-2], 3.9813262514556356)
        self.assertEqual(f.eta[-1], 0.0)
        
        
if __name__ == '__main__':
    unittest.main()