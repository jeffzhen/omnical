import omnical.info as Oi, omnical.calib as Oc, omnical._omnical as _O
#import omnical.calibration_omni as omni
import numpy as np, numpy.linalg as la
import os, unittest

redinfo_psa32 = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
#infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/redundantinfo_test.bin'
infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/calib_test_redinfo.npz'
testdata = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/calib_test_data_%02d.npz'

VERBOSE = False

class TestMethods(unittest.TestCase):
    def setUp(self):
        self.info = Oi.RedundantInfoLegacy(filename=redinfo_psa32, txtmode=True)
    def test_pack_calpar(self):
        calpar = np.zeros((2,3,Oc.calpar_size(self.info.nAntenna, len(self.info.ublcount))), dtype=np.float32)
        self.assertTrue(np.all(Oc.pack_calpar(self.info,calpar) == 0))
        self.assertRaises(AssertionError, Oc.pack_calpar, self.info, calpar[...,:-1])
        bp = np.array([[1+2j,3+4j,5+6j],[2+1j,4+3j,6+5j]])
        amp,phs = np.log10(np.abs(bp)), np.angle(bp)
        gains = {0:bp}
        Oc.pack_calpar(self.info,calpar,gains=gains)
        self.assertTrue(np.allclose(calpar[...,3+0], amp))
        self.assertTrue(np.allclose(calpar[...,32+3+0],phs))
        calpar *= 0
        gains = {1:bp[0]}
        Oc.pack_calpar(self.info,calpar,gains=gains)
        self.assertTrue(np.allclose(calpar[0,:,3+1], amp[0]))
        self.assertTrue(np.allclose(calpar[1,:,3+1], amp[0]))
        self.assertTrue(np.allclose(calpar[0,:,32+3+1],phs[0]))
        self.assertTrue(np.allclose(calpar[1,:,32+3+1],phs[0]))
        vis = {(0,16):bp}
        Oc.pack_calpar(self.info,calpar,vis=vis)
        self.assertTrue(np.allclose(calpar[...,3+2*32+2*12], bp.real))
        self.assertTrue(np.allclose(calpar[...,3+2*32+2*12+1], bp.imag))
    def test_unpack_calpar(self):
        calpar = np.zeros((2,3,Oc.calpar_size(self.info.nAntenna, len(self.info.ublcount))), dtype=np.float32)
        m,g,v = Oc.unpack_calpar(self.info,calpar)
        self.assertEqual(m['iter'].shape, (2,3))
        self.assertTrue(np.all(m['iter'] == 0))
        self.assertTrue(np.all(m['chisq'] == 0))
        self.assertEqual(len(g), 32)
        for i in xrange(32):
            self.assertTrue(np.all(g[i] == 1)) # 1 b/c 10**0 = 1
        self.assertEqual(len(v), len(self.info.ublcount))
        ubls = {}
        for i,j in v:
            n = self.info.bl1dmatrix[i,j]
            ubls[self.info.bltoubl[n]] = n
        for u in xrange(len(self.info.ublcount)):
            self.assertTrue(ubls.has_key(u))
    def test_redcal(self):
        #check that logcal give 0 chi2 for all 20 testinfos
        for index in xrange(20):
            arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(index+1)+'_array_info.txt'
            c = Oc.RedundantCalibrator(56)
            c.compute_redundantinfo(arrayinfopath, tol=.1)
            info = c.Info
            npz = np.load(testdata % index)
            bls = [tuple(bl) for bl in npz['bls']]
            dd = dict(zip(bls, npz['vis']))
            m,g,v = Oc.redcal(dd, info, removedegen=True,maxiter=50,stepsize=.2,computeUBLFit=True,conv=1e-5,uselogcal=True)
            calparpath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(index+1)+'_calpar.txt'
            with open(calparpath) as f:
                rawinfo = [[float(x) for x in line.split()] for line in f]
            temp = np.array(rawinfo[:-1])
            correctcalpar = (np.array(temp[:,0]) + 1.0j*np.array(temp[:,1]))
            i = g.keys()[0]
            scalar = correctcalpar[i].real / g[i].real
            for i in xrange(56):
                if not g.has_key(i): continue
                self.assertAlmostEqual(np.abs(correctcalpar[i] - g[i] * scalar), 0, 4)
    def test_redcal_xtalk(self):
        antpos = np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]])
        d = {(1,2): np.array([[1.]], dtype=np.complex64), (2,3): np.array([[1.+1j]], dtype=np.complex64)}
        x = {(1,2): np.array([[0.]], dtype=np.complex64), (2,3): np.array([[0.+1j]], dtype=np.complex64)}
        reds = [[(1,2),(2,3)]]
        info = Oi.RedundantInfo(); info.init_from_reds(reds, antpos)
        m,g,v = Oc.redcal(d, info, xtalk=x, uselogcal=False)
        self.assertEqual(g[1][0,0], 1.)
        self.assertEqual(g[2][0,0], 1.)
        self.assertEqual(g[3][0,0], 1.)
        #2D array testing
        d = {(1,2): np.array([[1.,2],[3.,4]],dtype=np.complex64), (2,3): np.array([[1.+1j,2+2j],[3+3j,4+4j]],dtype=np.complex64)}
        x = {(1,2): np.array([[0.,2],[2.,3]],dtype=np.complex64), (2,3): np.array([[0.+1j,1+2j],[2+3j,3+4j]],dtype=np.complex64)}
        m,g,v = Oc.redcal(d, info, xtalk=x, uselogcal=False)
        self.assertEqual(g[1][0,0], 1.)
        self.assertEqual(g[2][0,0], 1.)
        self.assertEqual(g[3][0,0], 1.)
        self.assertEqual(m['res'][(2,3)][0][0],0.)


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
            start_ubl = 3 + 2*info['nAntenna']
            end_ubl = start_ubl + 2*len(info.ublcount)
            ublfit = log[0,0,start_ubl:end_ubl:2]+1.0j*log[0,0,start_ubl+1:end_ubl+1:2]
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
        self.assertTrue((linlist < loglist).all())     #chick that chi2 of lincal is smaller than chi2 of logcal

    def test_loglincal_stdev1D(self):
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
        std = 0.1
        std_bl = np.random.uniform(.2, 1, size=data.shape[-1]).astype('float32')

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
            noise = (std_bl * (np.random.normal(scale=std, size=data.shape) + 1.0j*np.random.normal(scale=std, size = data.shape))).astype('complex64')
            ndata = data + noise
            calibrator.logcal(ndata, np.zeros_like(ndata), stdev=std_bl, verbose=VERBOSE)
            calibrator.lincal(ndata, np.zeros_like(ndata), stdev=std_bl, verbose=VERBOSE)

            linchi2 = (calibrator.rawCalpar[0,0,2]/(info['At'].shape[1] - info['At'].shape[0])/(2*(std)**2))**0.5
            logchi2 = (calibrator.rawCalpar[0,0,1]/(info['At'].shape[1] - info['At'].shape[0])/(2*(std)**2))**0.5
            linlist[i] = linchi2
            loglist[i] = logchi2
        # print np.mean(loglist), np.std(loglist)
        # print np.mean(linlist), np.std(linlist)
        self.assertTrue(abs(np.mean(linlist)-1.0) < 0.01)        #check that chi2 of lincal is close enough to 1
        self.assertTrue(np.mean(linlist) < np.mean(loglist))     #chick that chi2 of lincal is smaller than chi2 of logcal


    def test_lincal_stdev3D(self):
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
        std = 0.1
        std_bl = np.random.uniform(.2, 1, size=data.shape).astype('float32')

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
            noise = (std_bl * (np.random.normal(scale=std, size=data.shape) + 1.0j*np.random.normal(scale=std, size = data.shape))).astype('complex64')
            ndata = data + noise
            calibrator.logcal(ndata, np.zeros_like(ndata), stdev=std_bl[0,0], verbose=VERBOSE) #3D stdev not allowed for logcal
            calibrator.lincal(ndata, np.zeros_like(ndata), stdev=std_bl, verbose=VERBOSE)

            linchi2 = (calibrator.rawCalpar[0,0,2]/(info['At'].shape[1] - info['At'].shape[0])/(2*(std)**2))**0.5
            logchi2 = (calibrator.rawCalpar[0,0,1]/(info['At'].shape[1] - info['At'].shape[0])/(2*(std)**2))**0.5
            linlist[i] = linchi2
            loglist[i] = logchi2
        # print np.mean(loglist), np.std(loglist)
        # print np.mean(linlist), np.std(linlist)
        self.assertTrue(abs(np.mean(linlist)-1.0) < 0.01)        #check that chi2 of lincal is close enough to 1
        self.assertTrue(np.mean(linlist) < np.mean(loglist))     #chick that chi2 of lincal is smaller than chi2 of logcal

if __name__ == '__main__':
    unittest.main()
