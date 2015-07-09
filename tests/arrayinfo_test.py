import omnical.info as Oi, omnical.calib as Oc, omnical._omnical as _O, omnical.arrayinfo as Oa
#import omnical.calibration_omni as omni
import numpy as np, numpy.linalg as la
import os, unittest, aipy

redinfo_psa32 = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/redundantinfo_test.bin'

VERBOSE = False

class TestArrayInfo(unittest.TestCase):
    def test_compute_reds(self):
        diff = np.zeros(20)
        for index in range(20):
            fileindex = index+1
            arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_array_info.txt'
            ai = Oa.ArrayInfoLegacy(56)
            info = ai.compute_redundantinfo(arrayinfopath)
            ublgps = {}
            for bl in ai.totalVisibilityUBL:
                gp = ai.totalVisibilityUBL[bl]
                ublgps[gp] = ublgps.get(gp,[]) + [bl]
            ublreds = ublgps.values()
            ai.read_arrayinfo(arrayinfopath)
            reds = ai.compute_reds()
            # filter for "valid visibilities"
            reds = ai.filter_reds(reds, bls=ai.totalVisibilityId_dic.keys(), ex_ants=list(ai.badAntenna))
            for gp1 in reds: # make sure red groups are same as ubl groups (although conj of ublgps is wrong)
                if len(gp1) == 0: continue
                try: i = ai.totalVisibilityUBL[gp1[0]]
                except(KeyError): i = ai.totalVisibilityUBL[gp1[0][::-1]]
                gp2 = ublreds[i]
                for bl in gp1:
                    self.assertTrue(bl in gp2 or bl[::-1] in gp2)
                for bl in gp2: self.assertTrue(bl in gp1 or bl[::-1] in gp1)
    def test_filter_reds(self):
        reds = [[(1,2),(2,3),(3,4),(4,5)],[(1,3),(2,4),(3,5)],[(1,4),(2,5)]]
        self.assertEqual(Oa.filter_reds(reds,bls=[(1,2),(3,4),(4,5),(1,3),(3,5),(1,4),(2,5)]),
                [[(1,2),(3,4),(4,5)],[(1,3),(3,5)],[(1,4),(2,5)]])
        self.assertEqual(Oa.filter_reds(reds,bls=[(1,2),(3,4),(4,5)]),
                [[(1,2),(3,4),(4,5)]])
        self.assertEqual(Oa.filter_reds(reds,ex_bls=[(1,2),(3,4)]),
                [[(2,3),(4,5)],[(1,3),(2,4),(3,5)],[(1,4),(2,5)]])
        self.assertEqual(Oa.filter_reds(reds,ex_bls=[(2,5)]),
                [[(1,2),(2,3),(3,4),(4,5)],[(1,3),(2,4),(3,5)]])
        self.assertEqual(Oa.filter_reds(reds,ants=[1,2,3,4]),
                [[(1,2),(2,3),(3,4)],[(1,3),(2,4)]])
        self.assertEqual(Oa.filter_reds(reds,ex_ants=[5]),
                [[(1,2),(2,3),(3,4)],[(1,3),(2,4)]])
        self.assertEqual(Oa.filter_reds(reds,ubls=[(1,2),(1,3)]), reds[:2])
        self.assertEqual(Oa.filter_reds(reds,ubls=[(1,2)]), reds[:1])
        self.assertEqual(Oa.filter_reds(reds,ubls=[(1,3)]), reds[1:2])
        self.assertEqual(Oa.filter_reds(reds,ex_ubls=[(1,3)]), reds[0:1]+reds[2:3])
    def test_compute_reds_tol(self):
        antpos = np.zeros((65,3), dtype=np.float)
        antpos[3]  = np.array([ 3373240.12832542,   540919.61      ,  5674198.83539668])
        antpos[49] = np.array([ 3373240.64412165,   540904.7       ,  5674198.52887711])
        antpos[10] = np.array([ 3373240.90712842,   540889.64      ,  5674198.38421393])
        antpos[64] = np.array([ 3373240.51728633,   540874.74      ,  5674198.51119087])
        reds = Oa.filter_reds(Oa.compute_reds(antpos, tol=1.0), ants=(3,10,49,64))
        self.assertNotEqual(reds, [[(49, 64), (3, 10)], [(3, 49), (10, 64), (49, 10)]])
        reds = Oa.filter_reds(Oa.compute_reds(antpos, tol=1.5), ants=(3,10,49,64))
        self.assertEqual(reds, [[(49, 64), (3, 10)], [(3, 49), (10, 64), (49, 10)]])

if __name__ == '__main__':
    unittest.main()
