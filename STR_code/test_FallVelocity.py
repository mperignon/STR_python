import unittest
from FallVelocity import FallVelocity

class test_FallVelocity(unittest.TestCase):
    
    def test_Re(self):
        f = FallVelocity(D = 0.1, verbose = False, save_output = False)
        f.run()
        self.assertEqual(f.Reynolds_number, 4.023244958985222)
        
    def test_Rf(self):
        f = FallVelocity(D = 0.1, verbose = False, save_output = False)
        f.run()
        self.assertEqual(f.dimensionless_fall_velocity, 0.18603632526569483)
        
    def test_vs(self):
        f = FallVelocity(D = 0.1, verbose = False, save_output = False)
        f.run()
        self.assertEqual(f.settling_velocity, 0.007484697078133418)
        
if __name__ == '__main__':
    unittest.main()