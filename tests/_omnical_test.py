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

class TestRedInfo(unittest.TestCase):
    def test_getset(self):
        i = _O.RedundantInfo()
        i.nAntenna = 3
        self.assertEqual(i.nAntenna, 3)
        i.nUBL = 3
        self.assertEqual(i.nUBL, 3)
        i.nCross = 3
        self.assertEqual(i.nCross, 3)
    def test_getset_int1d(self):
        i = _O.RedundantInfo()
        ints = n.array([1,2,3], dtype=n.int32)
        for k in ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','ublcount']:
            i.__setattr__(k, ints)
            self.assertTrue(n.all(i.__getattribute__(k) == ints))
    def test_getset_int2d(self):
        i = _O.RedundantInfo()
        ints = n.array([[1,2,3],[4,5,6]], dtype=n.int32)
        for k in ['bl2d','bl1dmatrix','A','B','Atsparse']:
            i.__setattr__(k, ints)
            self.assertTrue(n.all(i.__getattribute__(k) == ints))
    def test_getset_int3d(self):
        i = _O.RedundantInfo()
        ints = n.array([[[1,2],[4,5]],[[5,6],[7,8]]], dtype=n.int32)
        for k in ['ublindex','Btsparse']:
            i.__setattr__(k, ints)
            self.assertTrue(n.all(i.__getattribute__(k) == ints))
    def test_getset_float2d(self):
        i = _O.RedundantInfo()
        floats = n.array([[1,2],[4,5]], dtype=n.float32)
        for k in ['antloc','ubl','degenM','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:
            i.__setattr__(k, floats)
            self.assertTrue(n.all(i.__getattribute__(k) == floats))

if __name__ == '__main__':
    unittest.main()
