import unittest, omnical._omnical as _O
import random
import numpy as np
import aipy as ap
import numpy.linalg as la
import commands, os, time, math, ephem, shutil
import omnical.calibration_omni as omni

print "#Omnical Version %s#"%omni.__version__

class TestTreasure(unittest.TestCase):
    def test_IO(self):
        nTime = 3
        nFrequency = 5
        shutil.rmtree(os.path.dirname(os.path.realpath(__file__)) + '/test.treasure', ignore_errors = True)
        shutil.rmtree(os.path.dirname(os.path.realpath(__file__)) + '/test2.treasure', ignore_errors = True)
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.add_coin(('xx', np.array([0,2,3])))
        treasure.add_coin(('xx', np.array([1,2,3])))
        self.assertEqual(treasure.coin_name(('xx', np.array([1,2,3]))), os.path.dirname(os.path.realpath(__file__)) + '/test.treasure//xx1.coin')

        treasure2 = treasure.duplicate_treasure(os.path.dirname(os.path.realpath(__file__)) + '/test2.treasure')
        treasure.burn()

        treasure2.add_coin(('xx', np.array([1,2,3])))
        treasure2.add_coin(('xx', np.array([1,2,4])))
        self.assertEqual(treasure2.coin_name(('xx', np.array([1,2,4]))), os.path.dirname(os.path.realpath(__file__)) + '/test2.treasure//xx2.coin')
        self.assertEqual(treasure2.coinShape, (nTime, nFrequency, 10))
        treasure2.burn()

    def test_math(self):
        nTime = 4
        nFrequency = 2
        shutil.rmtree(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', ignore_errors = True)
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.add_coin(('xx', np.array([0,2,3])))
        treasure.update_coin(('xx', np.array([0,2,3])), (treasure.lsts + treasure.lsts[1] * (nTime/2. + .5))%(2*np.pi), np.outer(np.arange(nTime), np.arange(nFrequency)), np.ones((nTime, nFrequency)))
        predict_result = np.outer(np.roll(np.append([0], (np.arange(nTime - 1) + np.arange(1, nTime)) / 2.), nTime/2, axis = 0), np.arange(nFrequency))
        #print (treasure.lsts + treasure.lsts[1] * (nTime/2. + .5))%(2*np.pi), np.outer(np.arange(nTime), np.arange(nFrequency))
        #print treasure.get_coin(('xx', np.array([0,2,3]))).mean
        #print predict_result
        #print predict_result - treasure.get_coin(('xx', np.array([0,2,3]))).mean
        np.testing.assert_almost_equal(predict_result, treasure.get_coin(('xx', np.array([0,2,3]))).mean, decimal = 14)

    def test_probability(self):
        nTime = 10
        nFrequency = 1
        shutil.rmtree(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', ignore_errors = True)
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.add_coin(('xx', np.array([0,2,3])))
        treasure.add_coin(('xx', np.array([1,2,3])))
        nupdate = 4
        update_lsts = np.append((treasure.lsts[-nupdate/2:]+np.pi/2/nTime), (treasure.lsts[:nupdate/2]+np.pi/2/nTime))
        nan_prob = .1
        trials = 10000
        for i in range(int(trials/(1-nan_prob))):
            #print i
            vis_re = (np.random.randn(nupdate) * (np.arange(nupdate) + 1) + range(nupdate)).reshape(nupdate, 1)
            vis_im = (np.random.randn(nupdate) * (np.arange(nupdate) + 1) + range(nupdate)).reshape(nupdate, 1)
            epsilons = (np.arange(nupdate, dtype='float') + 1).reshape(nupdate, 1)
            if random.random() < nan_prob:
                vis_re[:nupdate/2] = vis_re[:nupdate/2] + np.nan
            if random.random() < nan_prob:
                vis_re[-nupdate/2:] = vis_re[-nupdate/2:] + np.nan
            treasure.update_coin(('xx', np.array([1,2,3])), update_lsts, vis_re + 1.j * vis_im, epsilons**2)
        #print epsilons**2
        c = treasure.get_coin(('xx', np.array([1,2,3])))

        #print c.count, c.mean, c.weighted_mean
        #print c.variance_re, c.variance_im
        #print c.weighted_variance

        self.assertTrue(abs(c.count[1] - trials) < 3 * trials**.5)
        self.assertTrue(abs(c.count[-1] - trials) < 3 * trials**.5)

        sigma1 = (1/16. * epsilons[-2]**2 + 9/16. * epsilons[-1]**2)**.5
        sigma2 = (1/16. * epsilons[0]**2 + 9/16. * epsilons[1]**2)**.5
        for var in [c.weighted_variance, c.variance_re, c.variance_im]:
            weighted_sigma = (var * trials)**.5
            #print weighted_sigma, sigma1, sigma2
            self.assertTrue(abs(weighted_sigma[1] - sigma1)/sigma1 < 3 * trials**-.5)
            self.assertTrue(abs(weighted_sigma[-1] - sigma2)/sigma2 < 3 * trials**-.5)

        self.assertTrue(abs(c.mean[1] - 2.75-2.75j) < 1.414 * 3 * sigma1 * trials**-.5)
        self.assertTrue(abs(c.weighted_mean[1] - 2.75-2.75j) < 1.414 * 3 * sigma1 * trials**-.5)
        self.assertTrue(abs(c.mean[-1] - .75-.75j) < 1.414 * 3 * sigma2 * trials**-.5)
        self.assertTrue(abs(c.weighted_mean[-1] - .75-.75j) < 1.414 * 3 * sigma2 * trials**-.5)
        treasure.burn()

if __name__ == '__main__':
    unittest.main()
