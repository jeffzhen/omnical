#!/usr/bin/env python

import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import optparse, sys
FILENAME = "speedtest.py"


sides = [4,8,12,16]#,20,24]
noises = [0.0] + (10.**np.arange(-2, -.5, .5)).tolist()
times = np.zeros((len(sides), len(noises), 5))
phase_noise = .1

nt = 10
nf = 2
for trust_period in [1, 3, 10]:

	opname = "trust" + str(trust_period) + '_timing.txt'

	for nside, side in enumerate(sides):
		nant = side**2
		print "Starting %i antennas"%nant
		calibrator = omni.RedundantCalibrator(nant)

		print "Reading readundant info"
		calibrator.read_redundantinfo(os.path.dirname(os.path.realpath(__file__)) + '/../results/%i.bin'%nant, verbose = False)
		calibrator.removeDegeneracy = True
		calibrator.removeAdditive = False
		calibrator.keepData = False
		calibrator.keepCalpar = False
		calibrator.convergePercent = 1e-6
		calibrator.maxIteration = 200
		calibrator.stepSize = .3
		calibrator.computeUBLFit = False
		calibrator.nTime = nt
		calibrator.nFrequency = nf
		calibrator.trust_period = trust_period
		DoF = len(calibrator.Info.crossindex) - calibrator.Info.nUBL - calibrator.Info.nAntenna

		data = np.zeros((nt, nf, max(calibrator.Info.subsetbl) + 1), dtype='complex64')

		g_amp = np.random.randn(nant) * phase_noise
		gg_amp = g_amp[calibrator.Info.bl2d[calibrator.Info.crossindex, 0]] + g_amp[calibrator.Info.bl2d[calibrator.Info.crossindex, 1]]
		g_phase = np.random.randn(nant) * phase_noise
		gg_phase = -g_phase[calibrator.Info.bl2d[calibrator.Info.crossindex, 0]] + g_phase[calibrator.Info.bl2d[calibrator.Info.crossindex, 1]]
		for th1 in [.1, .4, .7, .8, 1.3, 1.4]:
			for ph1 in [.2, .8, 1.9, 2.8, 3.5, 4.5, 5]:
				k1 = np.array([np.sin(th1) * np.cos(ph1), np.sin(th1) * np.sin(ph1), np.cos(th1)])
				vis_phase = (2*np.pi*calibrator.Info.ubl.dot(k1))[calibrator.Info.bltoubl] * calibrator.Info.reversed

				data[:, :, calibrator.Info.subsetbl[calibrator.Info.crossindex]] += np.cos(th1) **2 * np.exp(1.j * vis_phase)[None,None,:] 
		data = data / np.mean(np.abs(data))
		model = np.copy(data[0,0])
		data[:, :, calibrator.Info.subsetbl[calibrator.Info.crossindex]] *= np.exp(gg_amp + 1.j * gg_phase)[None,None,:] 
		

		timer = omni.Timer()
		for nnoise, noise in enumerate(noises):
			ndata = data + (noise / 1.41421 * (np.random.randn(*data.shape) + np.random.randn(*data.shape) * 1.j)).astype('complex64')
			real_data = np.copy(ndata[0,0])
			#print ndata.shape, ndata.dtype
			timer.tick(noise, mute = True)
			calibrator.computeUBLFit = True
			calibrator.logcal(ndata, np.zeros_like(ndata), verbose=True)
			times[nside, nnoise, 0], _ = timer.tick(noise)
			calibrator.lincal(ndata, np.zeros_like(ndata), verbose=True)
			times[nside, nnoise, 1], _ = timer.tick(noise)
			times[nside, nnoise, 2] = np.average(calibrator.rawCalpar[:, :, 0])
			if noise != 0:
				times[nside, nnoise, 3] = np.average(calibrator.rawCalpar[:, :, 2])/noise**2./DoF
			else:
				times[nside, nnoise, 3] = np.average(calibrator.rawCalpar[:, :, 2])/DoF
			cal_data = calibrator.get_calibrated_data(ndata)[0,0]
			
			#sanity check: if input perfect calibration parameter, chi^2 should be ~1
			if noise != 0:
				times[nside, nnoise, 4] = np.sum(np.abs(ndata - data)**2)/ (nt*nf) / noise**2. / DoF
			else:
				times[nside, nnoise, 4] = np.sum(np.abs(ndata - data)**2)/ (nt*nf) / DoF
			#omni.omniview(np.array([model, real_data, cal_data]), calibrator.Info)
	oppath = os.path.dirname(os.path.realpath(__file__)) + '/../results/' + opname 
	while os.path.isfile(oppath):
		oppath += '_'
	np.savetxt(oppath, times.reshape((len(times) * len(times[0]), len(times[0,0]))))
