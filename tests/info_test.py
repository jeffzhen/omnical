import unittest, omnical.info as Oi, omnical._omnical as _O
import omnical.calibration_omni as omni
import numpy as np
import os

class TestRedundantInfo(unittest.TestCase):
    def test_init(self):
        i = Oi.RedundantInfo(threshold=64)
        self.assertEqual(i.threshold, 64)
        class RedInf(Oi.RedundantInfo):
            testcase = False
            def fromfile(self, filename, verbose=False, preview_only=False):
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
        self.assertTrue(i.compare(i, verbose=True))
    def test_fromfiletxt(self):
        filename = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
        i1 = omni.read_redundantinfo_txt(filename)
        i1 = omni.RedundantInfo(i1)
        i2 = Oi.RedundantInfo()
        i2.fromfile_txt(filename)
        self.assertTrue(i2.compare(i1,verbose=True))

if __name__ == '__main__':
    unittest.main()        
