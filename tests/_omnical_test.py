import unittest, omnical._omnical as _O
import numpy as np
import aipy as ap
import commands, os, time, math, ephem
import omnical.calibration_omni as omni

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
        self.assertTrue(np.all(calpar[:,:,:3+2*self.i.nAntenna] == 0)) # not great to be checking an initialization state
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna::2] == 0))
        self.assertTrue(np.all(calpar[:,:,3+2*self.i.nAntenna+1::2] == 0))
        self.assertTrue(np.all(additiveout == 0))

    def test_infoIO(self):
        correctinfo = omni.read_redundantinfo_txt(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt', verbose = False)
        infotestpath = './redundantinfo_test.bin'
        omni.write_redundantinfo(correctinfo, infotestpath, overwrite = True, verbose = False)
        Info = omni.RedundantInfo(infotestpath)
        self.assertTrue(omni.compare_info(correctinfo, Info.get_info(), tolerance = 1e-3))



    def test_all(self):
        ##FILENAME = "test.py"
        ######################################################################
        ##############Config parameters###################################
        ######################################################################
        ano = 'test'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
        uvfiles = [os.path.dirname(os.path.realpath(__file__)) + '/test.uv']
        wantpols = {'xx':-5}#, 'yy':-6}

        infopaths = {'xx':os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt', 'yy':os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'}
        arrayinfos = {'xx':os.path.dirname(os.path.realpath(__file__)) + '/../doc/arrayinfo_apprx_PAPER32.txt', 'yy':os.path.dirname(os.path.realpath(__file__)) + '/../doc/arrayinfo_apprx_PAPER32.txt'}

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
        calibrators = [omni.RedundantCalibrator(nant) for key in wantpols.keys()]
        for calibrator, key in zip(calibrators, wantpols.keys()):
            calibrator.compute_redundantinfo(arrayinfoPath = arrayinfos[key])
            #calibrator.write_redundantinfo(infoPath = './redundantinfo_test_' + key + '.txt', overwrite = True)
        ###start reading miriads################
        ##print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano
        data, t, timing, lst = omni.importuvs(uvfiles, calibrators[0].totalVisibilityId, wantpols, nTotalAntenna = 32,timingTolerance = 2*math.pi, init_mem = 5.e7)
        ##print FILENAME + " MSG:",  len(t), "slices read."

        ###raw calibration################
        if needrawcal:
            for p, key in zip(range(len(wantpols)), wantpols.keys()):
                rawcalpar = np.fromfile(rawpaths[key], dtype="complex64").reshape(nfreq, nant)
                data[p] = omni.apply_calpar(data[p], rawcalpar, calibrators[p].totalVisibilityId)

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
        for p, calibrator in zip(range(len(wantpols)), calibrators):
            calibrator.removeDegeneracy = removedegen
            calibrator.removeAdditive = removeadditive
            calibrator.keepData = keep_binary_data
            calibrator.keepCalpar = keep_binary_calpar
            calibrator.convergePercent = converge_percent
            calibrator.maxIteration = max_iter
            calibrator.stepSize = step_size
            calibrator.computeUBLFit = False


            calibrator.logcal(data[p], np.zeros_like(data[p]), verbose=True)
            calibrator.lincal(data[p], np.zeros_like(data[p]), verbose=True)

            calibrator.utctimes = timing
            calibrator.get_calibrated_data(data[p])
            calibrator.get_omnichisq()
            calibrator.get_omnigain()
            calibrator.get_omnifit()

        #########Test results############
        #calibrators[-1].rawCalpar.tofile(os.path.dirname(os.path.realpath(__file__)) + '/test.omnical')
        correctresult = np.fromfile(os.path.dirname(os.path.realpath(__file__)) + '/test.omnical', dtype = 'float32').reshape(14,203,165)#[:,:,3:]
        nanmask = ~np.isnan(np.sum(correctresult,axis=2))#mask the last dimension because when data contains some 0 and some -0, C++ code return various phasecalibration parameters on different systems, when all other numbers are nan. I do the summation to avoid it failing the euqality check when the input is trivially 0s.

        newresult = calibrators[-1].rawCalpar#[:,:,3:]

        #calibrators[-1].rawCalpar.tofile(os.path.dirname(os.path.realpath(__file__)) + '/test.omnical')
        np.testing.assert_almost_equal(correctresult[nanmask], newresult[nanmask])#decimal=

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
        for k in ['bl2d','bl1dmatrix','A','B']:#,'Atsparse']:
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
        for k in ['antloc','ubl','degenM','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:
            i.__setattr__(k, floats)
            self.assertTrue(np.all(i.__getattribute__(k) == floats))
    def test_readredundantinfo(self):
        i = _O.RedundantInfo()
        i.readredundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt')
        self.assertEqual(i.nAntenna, 32)
        self.assertEqual(i.subsetant.shape, (32,))


if __name__ == '__main__':
    unittest.main()
