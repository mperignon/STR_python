import unittest
from Backwater import Backwater

class test_Backwater(unittest.TestCase):
    
    def test_depth(self):
        f = Backwater(friction='Chezy',
                      qw = 5.7,
                      initdepth = 30.,
                      S = 0.00025,
                      total_x = 123000,
                      dx = 1000,
                      alpha = 8.1,
                      n = 2,
                      D90 = .91,
                      D50 = .27,
                      Cz = 22,
                      save_output = False)
        f.run()
        self.assertEqual(f.depth, 3.0136890806223846)
        
    def test_flow_velocity(self):
        f = Backwater(friction='Chezy',
                      qw = 5.7,
                      initdepth = 30.,
                      S = 0.00025,
                      total_x = 123000,
                      dx = 1000,
                      alpha = 8.1,
                      n = 2,
                      D90 = .91,
                      D50 = .27,
                      Cz = 22,
                      save_output = False)
        f.run()
        self.assertEqual(f.flow_velocity, 1.8913696295514468)
        
    def test_Fr(self):
        f = Backwater(friction='Chezy',
                      qw = 5.7,
                      initdepth = 30.,
                      S = 0.00025,
                      total_x = 123000,
                      dx = 1000,
                      alpha = 8.1,
                      n = 2,
                      D90 = .91,
                      D50 = .27,
                      Cz = 22,
                      save_output = False)
        f.run()
        self.assertEqual(f.Froude_number, 0.3478505426185218)
        
        
    def test_shear_stress(self):
        f = Backwater(friction='Chezy',
                      qw = 5.7,
                      initdepth = 30.,
                      S = 0.00025,
                      total_x = 123000,
                      dx = 1000,
                      alpha = 8.1,
                      n = 2,
                      D90 = .91,
                      D50 = .27,
                      Cz = 22,
                      save_output = False)
        f.run()
        self.assertEqual(f.shear_stress, 7.391072470226399)
        
if __name__ == '__main__':
    unittest.main()