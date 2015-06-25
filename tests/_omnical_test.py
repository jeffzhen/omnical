import unittest, omnical._omnical as _O
import omnical.info as Oi
import numpy as np
import numpy.linalg as la
import os, time, math

class TestMethods(unittest.TestCase):
    def setUp(self):
        #self.i = _O.RedundantInfo()
        self.i = Oi.RedundantInfoLegacy() # not ideal to require higher-level RedundantInfo object here, but no other way to read txt file
        self.i.fromfile_txt(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
    #    self.i.readredundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
    def test_phase(self):
        self.assertAlmostEqual(_O.phase(1.,0.), 0., 5)
        self.assertAlmostEqual(_O.phase(-1.,0.), np.pi, 5)

    def test_redcal_log(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.ublcount.size)),dtype='float32')
        additiveout = _O.redcal(data, calpar, self.i, additivein)
        self.assertTrue(np.all(calpar[:,:,:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 1))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_redcal_log_inplace(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        additiveout = np.zeros_like(data) + 1 + 1j
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.ublcount.size)),dtype='float32')
        _O.redcal(data, calpar, self.i, additivein, additiveout)
        self.assertTrue(np.all(calpar[:,:,:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 1))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_redcal_lin(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.ublcount.size)),dtype='float32')
        additiveout = _O.redcal(data, calpar, self.i, additivein, uselogcal=0)
        #print calpar[0,0,:3+2*self.i.nAntenna]
        self.assertTrue(np.all(calpar[:,:,:2] == 0))
        np.testing.assert_almost_equal(calpar[:,:,2], np.zeros((10,20)), 10)
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0)) # not great to be checking an initialization state
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_redcal_lin_inplace(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        additiveout = np.zeros_like(data) + 1 + 1j
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.ublcount.size)),dtype='float32')
        _O.redcal(data, calpar, self.i, additivein, additiveout, uselogcal=0)
        #print calpar[0,0,:3+2*self.i.nAntenna]
        self.assertTrue(np.all(calpar[:,:,:2] == 0))
        np.testing.assert_almost_equal(calpar[:,:,2], np.zeros((10,20)), 10)
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0)) # not great to be checking an initialization state
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_norm(self):
        d = np.zeros((2,3,4), dtype=np.float32)
        d[:,:,0] = 7
        #print _O.norm(d)
        np.testing.assert_array_equal(_O.norm(d), d.flatten()[:6])

class TestRedInfo(unittest.TestCase):
    def test_getset(self):
        i = _O.RedundantInfo()
        i.nAntenna = 3
        self.assertEqual(i.nAntenna, 3)
        ublcount = np.array([2,2,3], dtype=np.int32)
        i.ublcount = ublcount
        self.assertTrue(np.all(i.ublcount == ublcount))
        ublindex = np.array([15,16,4,3,2,1,0], dtype=np.int32)
        i.ublindex = ublindex
        self.assertTrue(np.all(i.ublindex == ublindex))
    def test_nosegfault(self):
        i = _O.RedundantInfo()
        for k in ['bl1dmatrix', 'bl2d', 'bltoubl', 'nAntenna', 'ublcount', 'ublindex']:
            i.__getattribute__(k)
    def test_getset_int1d(self):
        i = _O.RedundantInfo()
        ints = np.array([1,2,3], dtype=np.int32)
        for k in ['bltoubl','ublcount']:
            i.__setattr__(k, ints)
            self.assertTrue(np.all(i.__getattribute__(k) == ints))
    def test_getset_int2d(self):
        i = _O.RedundantInfo()
        ints = np.array([[1,2,3],[4,5,6]], dtype=np.int32)
        for k in ['bl2d','bl1dmatrix']:#,'A','B','Atsparse']:
            i.__setattr__(k, ints)
            self.assertTrue(np.all(i.__getattribute__(k) == ints))
    #def test_getset_int3d(self):
        #i = _O.RedundantInfo()
        #ints = np.array([[[1,2],[4,5]],[[5,6],[7,8]]], dtype=np.int32)
        #for k in ['ublindex','Btsparse']:
            #i.__setattr__(k, ints)
            #self.assertTrue(np.all(i.__getattribute__(k) == ints))
    def test_getset_float2d(self):
        i = _O.RedundantInfo()
        floats = np.array([[1,2],[4,5]], dtype=np.float32)
        for k in ['degenM','AtAi','BtBi']:#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:
            i.__setattr__(k, floats)
            self.assertTrue(np.all(i.__getattribute__(k) == floats))
    def test_readredundantinfo(self):
        i = Oi.RedundantInfoLegacy()
        #i = _O.RedundantInfo()
        i.fromfile_txt(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
        #i.readredundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
        self.assertEqual(i.nAntenna, 32)
        self.assertEqual(i.subsetant.shape, (32,))
    def test_getset_ublindex(self):
        ublcount = np.array([2,1], dtype=np.int32)
        ublindex = np.array([0, 2, 1], dtype=np.int32)
        i = _O.RedundantInfo()
        self.assertEqual(i.ublindex.size, 0)
        def f(val): i.ublindex = val
        self.assertRaises(ValueError, f, ublindex) # gotta set ublcount first
        i.ublcount = ublcount
        i.ublindex = ublindex
        self.assertTrue(np.all(i.ublindex == ublindex))
        self.assertRaises(ValueError, f, ublindex.reshape((1,3)))
        self.assertRaises(ValueError, f, ublindex.astype(np.float))
        self.assertRaises(ValueError, f, ublindex[:-1])
        self.assertRaises(ValueError, f, ublindex.reshape((3,1)))

if __name__ == '__main__':
    unittest.main()
