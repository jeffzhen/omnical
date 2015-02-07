import unittest, omnical._omnical as _O
import random
import numpy as np
import aipy as ap
import numpy.linalg as la
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
print "#Omnical Version %s#"%omni.__version__
class TestMethods(unittest.TestCase):
    def setUp(self):
        self.i = _O.RedundantInfo()
        self.i.readredundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')

    def test_phase(self):
        self.assertAlmostEqual(_O.phase(1.,0.), 0., 5)
        self.assertAlmostEqual(_O.phase(-1.,0.), np.pi, 5)

    def test_redcal_log(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.nUBL)),dtype='float32')
        additiveout = _O.redcal(data, calpar, self.i, additivein)
        self.assertTrue(np.all(calpar[:,:,:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 1))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_redcal_lin(self):
        data = np.ones((10,20,32*33/2), dtype=np.complex64)
        additivein = np.zeros_like(data)
        calpar = np.zeros((10,20,3+2*(self.i.nAntenna+self.i.nUBL)),dtype='float32')
        additiveout = _O.redcal(data, calpar, self.i, additivein, uselogcal=0)
        #print calpar[0,0,:3+2*self.i.nAntenna]
        self.assertTrue(np.all(calpar[:,:,:2] == 0))
        np.testing.assert_almost_equal(calpar[:,:,2], np.zeros((10,20)), 10)
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0))
        self.assertTrue(np.all(calpar[:,:,3:3+2*self.i.nAntenna] == 0)) # not great to be checking an initialization state
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_infoIO(self):
        correctinfo = omni.read_redundantinfo_txt(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt', verbose = False)
        infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/redundantinfo_test.bin'
        omni.write_redundantinfo(correctinfo, infotestpath, overwrite = True, verbose = False)
        Info = omni.RedundantInfo(infotestpath)
        self.assertTrue(omni.compare_info(correctinfo, Info.get_info(), tolerance = 1e-3))

    def test_testinfo_logcal(self):
        #check that logcal give 0 chi2 for all 20 testinfos
        diff = np.zeros(20)
        for index in range(20):
            fileindex = index+1   #filenames have index that start with 1
            ####Import arrayinfo##########################
            arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_array_info.txt'
            nant = 56
            calibrator = omni.RedundantCalibrator(nant)
            calibrator.compute_redundantinfo(arrayinfopath)
            info = calibrator.Info.get_info()
            ####Config parameters###################################
            removedegen = True
            removeadditive = False
            needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
            keep_binary_data = True
            keep_binary_calpar = True
            converge_percent = 0.00001
            max_iter = 50
            step_size = .2
            ####import data##################
            datapath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_data.txt'
            with open(datapath) as f:
                rawinfo = [[float(x) for x in line.split()] for line in f]
            data = np.array([i[0] + 1.0j*i[1] for i in rawinfo[:-1]],dtype = 'complex64')    #last element of rawinfo is empty
            data = data.reshape((1,1,len(data)))
            ####do calibration################
            calibrator.removeDegeneracy = removedegen
            calibrator.removeAdditive = removeadditive
            calibrator.keepData = keep_binary_data
            calibrator.keepCalpar = keep_binary_calpar
            calibrator.convergePercent = converge_percent
            calibrator.maxIteration = max_iter
            calibrator.stepSize = step_size
            calibrator.computeUBLFit = True

            calibrator.logcal(data, np.zeros_like(data), verbose=True)
            log = np.copy(calibrator.rawCalpar)
            log = np.copy(calibrator.rawCalpar)
            ampcal = log[0,0,3:info['nAntenna']+3]
            phasecal = log[0,0,info['nAntenna']+3: info['nAntenna']*2+3]
            calpar = 10**(ampcal)*np.exp(1.0j*phasecal)
            ublfit = log[0,0,3+2*info['nAntenna']::2]+1.0j*log[0,0,3+2*info['nAntenna']+1::2]
            ####import real calibration parameter
            calparpath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_calpar.txt'
            with open(calparpath) as f:
                rawinfo = [[float(x) for x in line.split()] for line in f]
            temp = np.array(rawinfo[:-1])
            correctcalpar = (np.array(temp[:,0]) + 1.0j*np.array(temp[:,1]))[info['subsetant']]
            ###compare calpar with correct calpar
            overallfactor = np.real(np.mean(ublfit))**0.5
            diffnorm = la.norm(calpar*overallfactor - correctcalpar)
            diff[index] = la.norm(diffnorm)
        self.assertAlmostEqual(la.norm(diff), 0, 4)

    def test_testinfo_lincal(self):
        fileindex = 3      #use the 3rd file to do the test, can also change this to any number from 1 to 20
        length = 100
        loglist = np.zeros(length)
        linlist = np.zeros(length)

        ####import arrayinfo################
        arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_array_info.txt'
        nant = 56
        calibrator = omni.RedundantCalibrator(nant)
        calibrator.compute_redundantinfo(arrayinfopath)
        info = calibrator.Info.get_info()

        ####Config parameters###################################
        removedegen = True
        removeadditive = False
        needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
        keep_binary_data = True
        keep_binary_calpar = True
        converge_percent = 0.00001
        max_iter = 50
        step_size = .2

        ####import data##################
        datapath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_data.txt'
        with open(datapath) as f:
            rawinfo = [[float(x) for x in line.split()] for line in f]
        data = np.array([i[0] + 1.0j*i[1] for i in rawinfo[:-1]],dtype = 'complex64')    #last element of rawinfo is empty
        truedata = data.reshape((1,1,len(data)))
        std = 0.1

        ####do calibration################
        calibrator.removeDegeneracy = removedegen
        calibrator.removeAdditive = removeadditive
        calibrator.keepData = keep_binary_data
        calibrator.keepCalpar = keep_binary_calpar
        calibrator.convergePercent = converge_percent
        calibrator.maxIteration = max_iter
        calibrator.stepSize = step_size
        calibrator.computeUBLFit = True

        for i in range(length):
            noise = (np.random.normal(scale = std, size = data.shape) + 1.0j*np.random.normal(scale = std, size = data.shape)).astype('complex64')
            data = truedata + noise
            calibrator.logcal(data, np.zeros_like(data), verbose=True)
            calibrator.lincal(data, np.zeros_like(data), verbose=True)

            linchi2 = (calibrator.rawCalpar[0,0,2]/(info['A'].shape[0] - info['A'].shape[1])/(2*std**2))**0.5
            logchi2 = (calibrator.rawCalpar[0,0,1]/(info['A'].shape[0] - info['A'].shape[1])/(2*std**2))**0.5
            linlist[i] = linchi2
            loglist[i] = logchi2
        self.assertTrue(abs(np.mean(linlist)-1.0) < 0.01)        #check that chi2 of lincal is close enough to 1
        self.assertTrue(np.mean(linlist) < np.mean(loglist))     #chick that chi2 of lincal is smaller than chi2 of logcal


    def test_norm(self):
        d = np.zeros((2,3,4), dtype=np.float32)
        d[:,:,0] = 7
        #print _O.norm(d)
        np.testing.assert_array_equal(_O.norm(d), d.flatten()[:6])

class TestUV(unittest.TestCase):
    def test_all(self):
        ##FILENAME = "test.py"
        ######################################################################
        ##############Config parameters###################################
        ######################################################################
        ano = 'test'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
        uvfiles = [os.path.dirname(os.path.realpath(__file__)) + '/test.uv']
        wantpols = {'xx':-5, 'yy':-6}

        #infopaths = {'xx':os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt', 'yy':os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'}
        arrayinfos = {'xx':os.path.dirname(os.path.realpath(__file__)) + '/../doc/arrayinfo_apprx_PAPER32_badUBLpair.txt', 'yy':os.path.dirname(os.path.realpath(__file__)) + '/../doc/arrayinfo_apprx_PAPER32_badUBLpair.txt'}

        oppath = os.path.dirname(os.path.realpath(__file__)) + '/results/'
        if not os.path.isdir(oppath):
            os.makedirs(oppath)
        removedegen = True
        removeadditive = False

        needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
        rawpaths = {'xx':os.path.dirname(os.path.realpath(__file__)) + "/testrawphasecalparrad_xx", 'yy':os.path.dirname(os.path.realpath(__file__)) + "/testrawphasecalparrad_yy"}



        keep_binary_data = True
        keep_binary_calpar = True


        converge_percent = 0.000001
        max_iter = 1000
        step_size = .3

        ######################################################################
        ######################################################################
        ######################################################################

        ########Massage user parameters###################################
        oppath += '/'


        ####get some info from the first uvfile   ################
        uv=ap.miriad.UV(uvfiles[0])
        nfreq = uv.nchan;
        nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
        startfreq = uv['sfreq']
        dfreq = uv['sdf']
        del(uv)

        ####create redundant calibrators################
        #calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
        calibrators = {}
        for key in wantpols.keys():
            calibrators[key] = omni.RedundantCalibrator(nant)
            calibrators[key].compute_redundantinfo(arrayinfoPath = arrayinfos[key])
            #calibrator.write_redundantinfo(infoPath = './redundantinfo_test_' + key + '.txt', overwrite = True)
        ###start reading miriads################
        ##print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano
        data, t, timing, lst, _ = omni.importuvs(uvfiles, wantpols, totalVisibilityId=calibrators[wantpols.keys()[0]].totalVisibilityId, nTotalAntenna = 32,timingTolerance = 2*math.pi, init_mem = 5.e7)
        ##print FILENAME + " MSG:",  len(t), "slices read."

        ###raw calibration################
        if needrawcal:
            for p, key in enumerate(wantpols.keys()):
                rawcalpar = np.fromfile(rawpaths[key], dtype="complex64").reshape(nfreq, nant)
                data[p] = omni.apply_calpar(data[p], rawcalpar, calibrators[key].totalVisibilityId)

        #####Save various files read################
        ###np.savetxt('miriadextract_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')
        ##f = open(oppath + 'miriadextract_' + ano + "_localtime.dat",'w')
        ##for time in timing:
            ##f.write("%s\n"%time)
        ##f.close()
        ##f = open(oppath + 'miriadextract_' + ano + "_lsthour.dat",'w')
        ##for l in lst:
            ##f.write("%s\n"%l)
        ##f.close()


        ####calibrate################
        ##print FILENAME + " MSG: starting calibration."
        for p, key in enumerate(wantpols.keys()):
            calibrators[key].removeDegeneracy = removedegen
            calibrators[key].removeAdditive = removeadditive
            calibrators[key].keepData = keep_binary_data
            calibrators[key].keepCalpar = keep_binary_calpar
            calibrators[key].convergePercent = converge_percent
            calibrators[key].maxIteration = max_iter
            calibrators[key].stepSize = step_size
            calibrators[key].computeUBLFit = False


            calibrators[key].logcal(data[p], np.zeros_like(data[p]), verbose=False)
            additiveout = calibrators[key].lincal(data[p], np.zeros_like(data[p]), verbose=False)

            calibrators[key].utctimes = timing
            calibrators[key].get_calibrated_data(data[p])
            calibrators[key].get_omnichisq()
            calibrators[key].get_omnigain()
            calibrators[key].get_omnifit()

        #########Test results############
        correctresult = np.fromfile(os.path.dirname(os.path.realpath(__file__)) + '/test.omnical', dtype = 'float32').reshape(14,203,165)[:,:,:]
        nanmask = ~np.isnan(np.sum(correctresult,axis=2))#mask the last dimension because when data contains some 0 and some -0, C++ code return various phasecalibration parameters on different systems, when all other numbers are nan. I do the summation to avoid it failing the euqality check when the input is trivially 0s.
        calibrators['xx'].rawCalpar.tofile(os.path.dirname(os.path.realpath(__file__)) + '/results/new_result.omnical')
        omni.write_redundantinfo(calibrators['xx'].Info.get_info(), os.path.dirname(os.path.realpath(__file__)) + '/results/new_info.bin', overwrite = True)
        newresult = calibrators['xx'].rawCalpar[:,:,:]
        np.testing.assert_almost_equal(correctresult[:,:,1:3][nanmask], newresult[:,:,1:3][nanmask], decimal = 5)
        np.testing.assert_almost_equal(np.sum(np.abs(data[wantpols.keys().index('xx')] - calibrators['xx'].get_modeled_data())[:,:,calibrators['xx'].Info.subsetbl[calibrators['xx'].Info.crossindex]]**2, axis=2)[nanmask], newresult[:,:,2][nanmask], decimal = 5)
        np.testing.assert_almost_equal(np.sum(np.abs(additiveout)[:,:,calibrators['xx'].Info.crossindex]**2, axis=2)[nanmask], newresult[:,:,2][nanmask], decimal = 5)
        np.testing.assert_almost_equal(correctresult[:,:,3:67][nanmask], newresult[:,:,3:67][nanmask], decimal = 5)
        np.testing.assert_almost_equal(np.sort(np.abs(correctresult[:,:,67:][nanmask])), np.sort(np.abs(newresult[:,:,67:][nanmask])), decimal = 5)

class TestTreasure(unittest.TestCase):
    def test_IO(self):
        nTime = 3
        nFrequency = 5
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.burn()
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.add_coin(np.array([0,2,3]))
        treasure.add_coin(np.array([1,2,3]))
        self.assertEqual(treasure.get_coin_index(np.array([1,2,3])), 1)

        treasure2 = treasure.duplicate_treasure(os.path.dirname(os.path.realpath(__file__)) + '/test2.treasure')
        treasure.burn()

        treasure2.add_coin(np.array([1,2,3]))
        treasure2.add_coin(np.array([1,2,4]))
        self.assertEqual(treasure2.get_coin_index(np.array([1,2,4])), 2)
        self.assertEqual(treasure2.coinShape, (nTime, nFrequency, 10))
        treasure2.burn()

    def test_math(self):
        nTime = 10
        nFrequency = 1
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.burn()
        treasure = omni.Treasure(os.path.dirname(os.path.realpath(__file__)) + '/test3.treasure', nlst = nTime, nfreq = nFrequency)
        treasure.add_coin(np.array([0,2,3]))
        treasure.add_coin(np.array([1,2,3]))
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
            treasure.update_coin(np.array([1,2,3]), update_lsts, vis_re + 1.j * vis_im, epsilons**2)
        #print epsilons**2
        c = treasure.get_coin(np.array([1,2,3]))

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
        ints = np.array([1,2,3], dtype=np.int32)
        for k in ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','ublcount']:
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
        for k in ['antloc','ubl','degenM','AtAi','BtBi']:#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:
            i.__setattr__(k, floats)
            self.assertTrue(np.all(i.__getattribute__(k) == floats))
    def test_readredundantinfo(self):
        i = _O.RedundantInfo()
        i.readredundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
        self.assertEqual(i.nAntenna, 32)
        self.assertEqual(i.subsetant.shape, (32,))


if __name__ == '__main__':
    unittest.main()
