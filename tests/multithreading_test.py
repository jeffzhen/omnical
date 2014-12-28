import unittest, omnical._omnical as _O
import random
import numpy as np
import aipy as ap
import numpy.linalg as la
import commands, os, time, math, ephem, multiprocessing, sys, copy
import omnical.calibration_omni as omni
import matplotlib.pyplot as plt
#def _f(m, m0, info, m1, ID, q):

	#time.sleep(info.nAntenna - 30)
	#m2 = np.copy(m)
	#m3 = 2*m2 + m2**2 + 3
	##print _O.phase(info.nAntenna,1)
	##_O.redcal(m, m0, info, m1)
	#print "putting onto ", ID,
	#q.put((ID,m[0]), block=False)
	#print "Done", ID
	#return 0


#timer = omni.Timer()
#ts={}
#matrix = np.zeros((2000,200,600), dtype='complex64')
#calpar = np.zeros((2000,200,600), dtype='float32')
#nthread = 10

#calibrator = omni.RedundantCalibrator(32)
#calibrator.compute_redundantinfo(arrayinfoPath = os.path.dirname(os.path.realpath(__file__)) + '/../doc/arrayinfo_apprx_PAPER32_badUBLpair.txt')
#q = multiprocessing.Queue()
#for i in range(nthread):
	#info = omni.RedundantInfo(calibrator.Info.get_info())
	#ts[i] = multiprocessing.Process(target = _f, args = (matrix[:, i::nthread, calibrator.Info.subsetbl], calpar[:, i::nthread, :3 + 2*(calibrator.Info.nAntenna + calibrator.Info.nUBL)], info, matrix[:, i::nthread, calibrator.Info.subsetbl], i, q))
	##ts[i] = threading.Thread(target = _O.redcal, args = (matrix[:, i::nthread, calibrator.Info.subsetbl], calpar[:, i::nthread, :3 + 2*(calibrator.Info.nAntenna + calibrator.Info.nUBL)], calibrator.Info, matrix[:, i::nthread, calibrator.Info.subsetbl]))
#for i in range(nthread):
	#print "starting", i
	#ts[i].start()
#for i in range(nthread):
	#print "collecting", i
	#ts[i].join()
	#print "joined", i
	##op = q.get()
	##print "Got from queue", op[0], op[1]
#timer.tick()

#exit()

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
	data = np.concatenate([data[p] for i in range(1)], axis = 0)
	calibrator.removeDegeneracy = removedegen
	calibrator.removeAdditive = removeadditive
	calibrator.keepData = keep_binary_data
	calibrator.keepCalpar = keep_binary_calpar
	calibrator.convergePercent = converge_percent
	calibrator.maxIteration = max_iter
	calibrator.stepSize = step_size
	calibrator.computeUBLFit = False
	#print data[0,120:132,0]
	times = []
	for nthread in range(1, 26):
		timer = omni.Timer()
		calibrator.logcal(data, np.zeros_like(data), nthread = nthread, verbose=False)
		#timer.tick()
		calibrator.lincal(data, np.zeros_like(data), nthread = nthread, verbose=False)
		times += [list((nthread,) + timer.tick(nthread))]
times = np.array(times)
plt.plot(times[:, 1])
plt.plot(times[1, 1] * 2/np.array(range(1,26)))
plt.show()
