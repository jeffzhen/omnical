#!/usr/bin/env python

import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
FILENAME = "speedtest.py"

sides = [4,8,12,16,20]#,24]
noises = 10.**np.arange(-2, -.5, .5)
times = np.zeros((len(sides), len(noises), 2))

nt = 10
nf = 2

for nside, side in enumerate(sides):
    nant = side**2
    print "Starting %i antennas"%nant
    calibrator = omni.RedundantCalibrator(nant)

    print "Reading readundant info"
    calibrator.read_redundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../results/%i.bin'%nant, verbose = False)
    calibrator.removeDegeneracy = False
    calibrator.removeAdditive = False
    calibrator.keepData = False
    calibrator.keepCalpar = False
    calibrator.convergePercent = 1e-6
    calibrator.maxIteration = 200
    calibrator.stepSize = .3
    calibrator.computeUBLFit = False
    calibrator.nTime = nt
    calibrator.nFrequency = nf

    data = np.zeros((nt, nf, max(calibrator.Info.subsetbl) + 1), dtype='complex64')
    th1, ph1 = .5, .5
    k1 = np.array([np.sin(th1) * np.cos(ph1), np.sin(th1) * np.sin(ph1), np.cos(th1)])
    data[:, :, calibrator.Info.subsetbl[calibrator.Info.crossindex]] = np.exp((2*np.pi*1.j*calibrator.Info.ubl.dot(k1))[calibrator.Info.bltoubl] * calibrator.Info.reversed)[None,None,:] 
    #omni.omniview(data[0,0], calibrator.Info)

    timer = omni.Timer()
    for nnoise, noise in enumerate(noises):
        ndata = data + (noise / 1.414 * (np.random.randn(*data.shape) + np.random.randn(*data.shape) * 1.j)).astype('complex64')
        #print ndata.shape, ndata.dtype
        timer.tick(noise, mute = True)
        calibrator.logcal(ndata, np.zeros_like(ndata), verbose=True)
        times[nside, nnoise, 0], _ = timer.tick(noise)
        calibrator.lincal(ndata, np.zeros_like(ndata), verbose=True)
        times[nside, nnoise, 1], _ = timer.tick(noise)

np.savetxt(os.path.dirname(os.path.realpath(__file__)) + '/../results/timing.txt', times.reshape((len(times) * len(times[0]), len(times[0,0]))))
