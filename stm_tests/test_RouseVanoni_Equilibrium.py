import unittest
from RouseVanoni_Equilibrium import RouseVanoni_Equilibrium

class test_RouseVanoni_Equilibrium(unittest.TestCase):
    
    def test_profile1(self):
        f = RouseVanoni_Equilibrium(b = 0.05, u_star = 0.2, vs = 3., save_output = False)
        self.assertEqual(f.C[1], 0.75562851889846694)
        
    def test_profile2(self):
        f = RouseVanoni_Equilibrium(b = 0.05, u_star = 0.2, vs = 3., save_output = False)
        self.assertEqual(f.C[-2], 0.045540468694966078)
        
    def test_profile3(self):
        f = RouseVanoni_Equilibrium(b = 0.05, u_star = 0.2, vs = 3., save_output = False)
        self.assertEqual(f.C[10], 0.30745820244404254)
        
if __name__ == '__main__':
    unittest.main()