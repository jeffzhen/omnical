import omnical.info as Oi, omnical.calib as Oc, omnical._omnical as _O
#import omnical.calibration_omni as omni
import numpy as np, numpy.linalg as la
import os, unittest

redinfo_psa32 = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
#infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/redundantinfo_test.bin'
infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/calib_test_redinfo.npz'
testdata = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/calib_test_data_%02d.npz'

VERBOSE = False

class TestRedCal(unittest.TestCase):
    #def setUp(self):
    #    self.i = Oi.RedundantInfo()
    #    self.i.fromfile_txt(redinfo_psa32)
    def tearDown(self):
        if os.path.exists(infotestpath): os.remove(infotestpath)

    def test_large_info_IO(self):
        calibrator = Oc.RedundantCalibrator(150)
        calibrator.compute_redundantinfo()
        calibrator.write_redundantinfo(infotestpath, verbose=VERBOSE)
        info2 = Oi.RedundantInfo(filename=infotestpath)
        self.assertEqual(calibrator.Info.nAntenna, info2.nAntenna)
        self.assertEqual(calibrator.Info.nBaseline, info2.nBaseline)
        self.assertEqual(calibrator.Info.get_reds(), info2.get_reds())
        os.remove(infotestpath)

    def test_logcal(self):
        #check that logcal give 0 chi2 for all 20 testinfos
        diff = np.zeros(20)
        for index in range(20):
            arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(index+1)+'_array_info.txt'
            calibrator = Oc.RedundantCalibrator(56)
            calibrator.compute_redundantinfo(arrayinfopath, tol=.1)
            if False: # XXX this was to migrate files so they include bl order w/ data
                _info = calibrator.Info # XXX needs to have been initialized the old way (w/o reds)
                datapath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(index+1)+'_data.txt'
                with open(datapath) as f:
                    rawinfo = [[float(x) for x in line.split()] for line in f]
                data = np.array([i[0] + 1.0j*i[1] for i in rawinfo[:-1]],dtype = 'complex64') #last element is empty
                data = data.reshape((1,1,len(data)))
                dd = _info.make_dd(data)
                np.savez('calib_test_data_%02d.npz' % index, bls=np.array(dd.keys()), vis=np.array(dd.values()))
            info = calibrator.Info
            npz = np.load(testdata % index)
            bls = [tuple(bl) for bl in npz['bls']]
            dd = dict(zip(bls, npz['vis']))
            data = info.order_data(dd)
            ####do calibration################
            calibrator.removeDegeneracy = True
            calibrator.removeAdditive = False
            calibrator.keepData = True
            calibrator.keepCalpar = True
            calibrator.convergePercent = 1e-5
            calibrator.maxIteration = 50
            calibrator.stepSize = .2
            calibrator.computeUBLFit = True

            calibrator.logcal(data, np.zeros_like(data), verbose=VERBOSE)
            log = np.copy(calibrator.rawCalpar)
            ampcal = log[0,0,3:info['nAntenna']+3]
            phasecal = log[0,0,info['nAntenna']+3: info['nAntenna']*2+3]
            calpar = 10**(ampcal)*np.exp(1.0j*phasecal)
            ublfit = log[0,0,3+2*info['nAntenna']::2]+1.0j*log[0,0,3+2*info['nAntenna']+1::2]
            ####import real calibration parameter
            calparpath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(index+1)+'_calpar.txt'
            with open(calparpath) as f:
                rawinfo = [[float(x) for x in line.split()] for line in f]
            temp = np.array(rawinfo[:-1])
            correctcalpar = (np.array(temp[:,0]) + 1.0j*np.array(temp[:,1]))[info['subsetant']]
            ###compare calpar with correct calpar
            overallfactor = np.real(np.mean(ublfit))**0.5
            diffnorm = la.norm(calpar*overallfactor - correctcalpar)
            self.assertAlmostEqual(diffnorm, 0, 4)

    def test_lincal(self):
        fileindex = 3      #use the 3rd file to do the test, can also change this to any number from 1 to 20
        length = 100
        loglist = np.zeros(length)
        linlist = np.zeros(length)

        ####import arrayinfo################
        arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_array_info.txt'
        nant = 56
        calibrator = Oc.RedundantCalibrator(nant)
        calibrator.compute_redundantinfo(arrayinfopath)
        info = calibrator.Info
        npz = np.load(testdata % (fileindex-1))
        bls = [tuple(bl) for bl in npz['bls']]
        dd = dict(zip(bls, npz['vis']))
        data = info.order_data(dd)

        ####Config parameters###################################
        needrawcal = True #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
        std = 0.1

        ####do calibration################
        calibrator.removeDegeneracy = True
        calibrator.removeAdditive = False
        calibrator.keepData = True
        calibrator.keepCalpar = True
        calibrator.convergePercent = 1e-5
        calibrator.maxIteration = 50
        calibrator.stepSize = .2
        calibrator.computeUBLFit = True

        for i in range(length):
            noise = (np.random.normal(scale = std, size = data.shape) + 1.0j*np.random.normal(scale = std, size = data.shape)).astype('complex64')
            ndata = data + noise
            calibrator.logcal(ndata, np.zeros_like(ndata), verbose=VERBOSE)
            calibrator.lincal(ndata, np.zeros_like(ndata), verbose=VERBOSE)

            linchi2 = (calibrator.rawCalpar[0,0,2]/(info['At'].shape[1] - info['At'].shape[0])/(2*std**2))**0.5
            logchi2 = (calibrator.rawCalpar[0,0,1]/(info['At'].shape[1] - info['At'].shape[0])/(2*std**2))**0.5
            linlist[i] = linchi2
            loglist[i] = logchi2
        self.assertTrue(abs(np.mean(linlist)-1.0) < 0.01)        #check that chi2 of lincal is close enough to 1
        self.assertTrue(np.mean(linlist) < np.mean(loglist))     #chick that chi2 of lincal is smaller than chi2 of logcal

if __name__ == '__main__':
    unittest.main()
