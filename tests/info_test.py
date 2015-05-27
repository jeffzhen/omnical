import unittest, omnical.info as Oi, omnical._omnical as _O
import omnical.calibration_omni as omni
import numpy as np
import os

redinfo_psa32 = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
tmpfile = 'aiwvlkasfdlk'
tmpnpz = 'xcmowtpqpow.npz'
VERBOSE = False

class TestRedInfo(unittest.TestCase):
    def tearDown(self):
        if os.path.exists(tmpfile): os.remove(tmpfile)
        if os.path.exists(tmpnpz): os.remove(tmpnpz)
    def test_init(self):
        i = Oi.RedundantInfo(threshold=64)
        self.assertEqual(i.threshold, 64)
        class RedInf(Oi.RedundantInfo):
            testcase = False
            def fromfile(self, filename, verbose=VERBOSE, preview_only=False):
                self.testcase = True
        i = RedInf(filename='whatever')
        self.assertTrue(i.testcase)
    def test_basic_getsetitem(self):
        i = Oi.RedundantInfo()
        i.nAntenna = 32
        self.assertEqual(i.nAntenna, 32)
        self.assertEqual(i['nAntenna'], 32)
        i['nAntenna'] = 64
        self.assertEqual(i.nAntenna, 64)
        self.assertEqual(i['nAntenna'], 64)
    def test_compare(self):
        i = Oi.RedundantInfo()
        self.assertTrue(i.compare(i, verbose=VERBOSE))
    def test_fromfiletxt(self):
        i1 = omni.read_redundantinfo_txt(redinfo_psa32)
        i1 = omni.RedundantInfo(i1)
        i2 = Oi.RedundantInfo()
        i2.fromfile_txt(redinfo_psa32)
        self.assertTrue(i2.compare(i1,verbose=VERBOSE))
    def test_tofromarray(self):
        i1 = Oi.RedundantInfo()
        i1.fromfile_txt(redinfo_psa32)
        i2 = Oi.RedundantInfo()
        d = i1.to_array(verbose=VERBOSE)
        # XXX yuck
        d = np.array(d)
        markerindex = np.where(d == Oi.MARKER)[0]
        d = np.array([np.array(d[markerindex[i]+1:markerindex[i+1]]) for i in range(len(markerindex)-1)])
        i2.from_array(d)
        i2.update()
        self.assertTrue(i1.compare(i2, verbose=VERBOSE))
    def test_tofromfile(self):
        i1 = Oi.RedundantInfo()
        i1.fromfile_txt(redinfo_psa32)
        i2 = Oi.RedundantInfo()
        i1.tofile(tmpfile)
        i2.fromfile(tmpfile)
        self.assertTrue(i1.compare(i2, verbose=VERBOSE))
    def test_tofromnpz(self):
        i1 = Oi.RedundantInfo()
        i1.fromfile_txt(redinfo_psa32)
        i2 = Oi.RedundantInfo()
        i1.to_npz(tmpnpz)
        i2.from_npz(tmpnpz)
        self.assertTrue(i1.compare(i2, verbose=VERBOSE))

if __name__ == '__main__':
    unittest.main()        
