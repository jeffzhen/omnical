import unittest, omnical._omnical as _O
import numpy as n

class TestImport(unittest.TestCase):
    def test_phase(self):
        self.assertAlmostEqual(_O.phase(1.,0.), 0., 5)
        self.assertAlmostEqual(_O.phase(-1.,0.), n.pi, 5)
    def test_norm(self):
        d = n.zeros((20,), dtype=n.float32)
        d[0] = 7
        self.assertEqual(_O.norm(d), d[0])

if __name__ == '__main__':
    unittest.main()
