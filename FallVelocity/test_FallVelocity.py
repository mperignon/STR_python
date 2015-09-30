import unittest
from FallVelocity import FallVelocity

class test_FallVelocity(unittest.TestCase):
    
    def test_Re(self):
        f = FallVelocity(0.1, verbose = False)
        self.assertEqual(f.Re, 4.023244958985222)
        
    def test_Rf(self):
        f = FallVelocity(0.1, verbose = False)
        self.assertEqual(f.Rf, 0.1891915672463031)
        
    def test_vs(self):
        f = FallVelocity(0.1, verbose = False)
        self.assertEqual(f.vs, 0.007611640192062025)
        
if __name__ == '__main__':
    unittest.main()