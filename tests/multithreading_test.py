import unittest, omnical._omnical as _O
import random
import numpy as np
import aipy as ap
import numpy.linalg as la
import commands, os, time, math, ephem
import omnical.calibration_omni as omni

ano = 'test'##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
uvfiles = [os.path.dirname(os.path.realpath(__file__)) + '/test.uv']
wantpols = {'xx':-5}#, 'yy':-6}

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
	data = np.concatenate(list([data[p] for i in range(5)]), axis = 0)
	calibrator.removeDegeneracy = removedegen
	calibrator.removeAdditive = removeadditive
	calibrator.keepData = keep_binary_data
	calibrator.keepCalpar = keep_binary_calpar
	calibrator.convergePercent = converge_percent
	calibrator.maxIteration = max_iter
	calibrator.stepSize = step_size
	calibrator.computeUBLFit = False

	timer = omni.Timer()
	calibrator.logcal(data, np.zeros_like(data), nthread = None, verbose=True)
	timer.tick()
