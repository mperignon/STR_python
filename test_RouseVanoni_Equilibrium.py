import unittest
from RouseVanoni_Equilibrium import RouseVanoni_Equilibrium

class test_RouseVanoni_Equilibrium(unittest.TestCase):
    
    def test_profile(self):
        f = RouseVanoni_Equilibrium(b = 0.05,
                                    u_star = 0.2,
                                    settling_velocity = 0.03,
                                    save_output = False)
        f.run()
                                    
        self.assertEqual(f.RV_profile[1,1], 0.75562851889846694)
        self.assertEqual(f.RV_profile[-2,1], 0.04554046869496612)
        self.assertEqual(f.RV_profile[10,1], 0.30745820244404265)
        
if __name__ == '__main__':
    unittest.main()