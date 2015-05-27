import unittest, omnical.info as Oi, omnical.calib as Oc, omnical._omnical as _O
import omnical.calibration_omni as omni
import random
import numpy as np
import aipy as ap
import numpy.linalg as la
import commands, os, time, math, ephem, shutil

redinfo_psa32 = os.path.dirname(os.path.realpath(__file__)) + '/../doc/redundantinfo_PSA32.txt'
infotestpath = os.path.dirname(os.path.realpath(__file__)) + '/redundantinfo_test.bin'

VERBOSE = False

class TestRedundantCalibrator(unittest.TestCase):
    def setUp(self):
        self.i = Oi.RedundantInfo()
        self.i.fromfile_txt(redinfo_psa32)
    def tearDown(self):
        if os.path.exists(infotestpath): os.remove(infotestpath)

    def test_large_info_IO(self):
        calibrator = Oc.RedundantCalibrator(150)
        calibrator.compute_redundantinfo()
        calibrator.write_redundantinfo(infotestpath, verbose=VERBOSE)
        info2 = Oi.RedundantInfo(filename=infotestpath)
        self.assertTrue(calibrator.Info.compare(info2, tol=1e-3))
        os.remove(infotestpath)

    def test_testinfo_logcal(self):
        #check that logcal give 0 chi2 for all 20 testinfos
        diff = np.zeros(20)
        for index in range(20):
            fileindex = index+1   #filenames have index that start with 1
            ####Import arrayinfo##########################
            arrayinfopath = os.path.dirname(os.path.realpath(__file__)) + '/testinfo/test'+str(fileindex)+'_array_info.txt'
            nant = 56
            #calibrator = omni.RedundantCalibrator(nant)
            calibrator = Oc.RedundantCalibrator(nant)
            calibrator.compute_redundantinfo(arrayinfopath)
            #self.assertTrue(calibrator2.Info.compare(calibrator.Info, tol=1e-3))
            #info = calibrator.Info.get_info()
            info = calibrator.Info
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

            calibrator.logcal(data, np.zeros_like(data), verbose=VERBOSE)
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
        #calibrator = omni.RedundantCalibrator(nant)
        calibrator = Oc.RedundantCalibrator(nant)
        calibrator.compute_redundantinfo(arrayinfopath)
        #info = calibrator.Info.get_info()
        info = calibrator.Info

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
            calibrator.logcal(data, np.zeros_like(data), verbose=VERBOSE)
            calibrator.lincal(data, np.zeros_like(data), verbose=VERBOSE)

            linchi2 = (calibrator.rawCalpar[0,0,2]/(info['At'].shape[1] - info['At'].shape[0])/(2*std**2))**0.5
            logchi2 = (calibrator.rawCalpar[0,0,1]/(info['At'].shape[1] - info['At'].shape[0])/(2*std**2))**0.5
            linlist[i] = linchi2
            loglist[i] = logchi2
        self.assertTrue(abs(np.mean(linlist)-1.0) < 0.01)        #check that chi2 of lincal is close enough to 1
        self.assertTrue(np.mean(linlist) < np.mean(loglist))     #chick that chi2 of lincal is smaller than chi2 of logcal

if __name__ == '__main__':
    unittest.main()
