import unittest
from GSD_calculator import GSD_calculator

class test_GSD_calculator(unittest.TestCase):
    
    def test_grain_processing(self):
    
        datapts = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), 
               (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
    
        f = GSD_calculator(datapts = datapts)
        f.run()
        
        self.assertEqual(f.datapts[0][0], 0.0097643881366852696)
        
    def test_geometric_mean(self):
    
        datapts = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), 
               (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
    
        
        f = GSD_calculator(datapts = datapts)
        f.run()
        
        self.assertEqual(f.geometric_mean, 0.2729234998323391)
        
    def test_std(self):
    
        datapts = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), 
               (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
    
        
        f = GSD_calculator(datapts = datapts)
        f.run()
        
        self.assertEqual(f.standard_deviation, 2.1681747768957598)
        
        
    def test_grain_statistics(self):
    
        datapts = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), 
               (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
    
        
        f = GSD_calculator(datapts = datapts)
        f.run()
        grain_stats = f.characteristic_size([70,50,30])
        
        self.assertEqual(grain_stats[0], 0.39951691807960471)
        self.assertEqual(grain_stats[1], 0.28583147875448067)
        self.assertEqual(grain_stats[2], 0.19277635317599259)
        
if __name__ == '__main__':
    unittest.main()