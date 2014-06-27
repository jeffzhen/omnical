#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import calibration_omni as omni
import optparse, sys
FILENAME = "omnical_PSA64.py"

##########################Sub-class#############################
class RedundantCalibrator_PAPER(omni.RedundantCalibrator):
	def __init__(self, aa):
		nTotalAnt = len(aa)
		omni.RedundantCalibrator.__init__(self, nTotalAnt)
		self.aa = aa

	def compute_redundantinfo(self, badAntenna = [], badUBL = [], antennaLocationTolerance = 1e-6):
		self.antennaLocationTolerance = antennaLocationTolerance
		self.badAntenna = badAntenna
		self.badUBL = badUBL
		self.antennaLocation = np.zeros((self.nTotalAnt,3))
		for i in range(len(self.aa.ant_layout)):
			for j in range(len(self.aa.ant_layout[0])):
				self.antennaLocation[self.aa.ant_layout[i][j]] = np.array([i, j, 0])
		self.preciseAntennaLocation = np.array([ant.pos for ant in self.aa])
		omni.RedundantCalibrator.compute_redundantinfo(self)





######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

ap.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--tag', action = 'store', default = 'PSA64', help = 'tag name of this calculation')
o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
o.add_option('--skip', action = 'store_true', help = 'whether to skip data importing')
opts,args = o.parse_args(sys.argv[1:])
skip = opts.skip

ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
uvfiles = args
for uvf in uvfiles:
	if not os.path.isdir(uvf):
		uvfiles.remove(uvf)
if len(uvfiles) == 0:
	print "ERROR: No valid uv files detected in input. Exiting!"
	quit()

wantpols = {}
for p in opts.pol.split(','): wantpols[p] = ap.miriad.str2pol[p]
#wantpols = {'xx':ap.miriad.str2pol['xx']}#, 'yy':-6}#todo:


aa = ap.cal.get_aa(opts.cal, np.array([.15]))

badAntenna = [37]
badUBL = []

oppath = './results/'

infopaths = {'xx':oppath + 'redundantinfo_PSA64.txt', 'yy':oppath + 'redundantinfo_PSA64.txt'}
#arrayinfos = {'xx':'./arrayinfo_apprx_PAPER32.txt', 'yy':'./arrayinfo_apprx_PAPER32.txt'}


removedegen = True
if opts.add and opts.nadd > 0:
	removeadditive = True
	removeadditiveperiod = opts.nadd
else:
	removeadditive = False
	removeadditiveperiod = -1

needrawcal = False #if true, (generally true for raw data) you need to take care of having raw calibration parameters in float32 binary format freq x nant
#rawpaths = {'xx':"testrawphasecalparrad_xx", 'yy':"testrawphasecalparrad_yy"}



keep_binary_data = True
keep_binary_calpar = True


converge_percent = 0.001
max_iter = 20
step_size = .3

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
oppath += '/'
utcPath = oppath + 'miriadextract_' + ano + "_localtime.dat"
lstPath = oppath + 'miriadextract_' + ano + "_lsthour.dat"

####get some info from the first uvfile   ################
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
#nant = uv['nants']
#startfreq = uv['sfreq']
#dfreq = uv['sdf']
del(uv)


####create redundant calibrators################
#calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
calibrators = [RedundantCalibrator_PAPER(aa) for key in wantpols.keys()]
for calibrator, key in zip(calibrators, wantpols.keys()):
	#calibrator.compute_redundantinfo(badAntenna = badAntenna, badUBL = badUBL, antennaLocationTolerance = 1)
	#calibrator.write_redundantinfo(infoPath = oppath + 'redundantinfo_' + ano + '_' + key + '.txt', overwrite = True)
	calibrator.read_redundantinfo(infopaths[p])
	calibrator.dataPath = oppath + 'data_' + ano + '_' + key
	calibrator.tmpDataPath = calibrator.dataPath
	if removeadditive:
		calibrator.calparPath = oppath + 'data_' + ano + '_add' + str(removeadditiveperiod) + '_' + key + '.omnical'
	else:
		calibrator.calparPath = oppath + 'data_' + ano + '_' + key + '.omnical'


###start reading miriads################
if skip:
	with open(utcPath) as f:
		timing = f.readlines()
else:
	print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano
	data, t, timing, lst = omni.importuvs(uvfiles, calibrators[0].totalVisibilityId, wantpols)#, nTotalAntenna = len(aa))
	print FILENAME + " MSG:",  len(t), "slices read."
	f = open(utcPath,'w')
	for time in timing:
		f.write("%s\n"%time)
	f.close()
	f = open(lstPath,'w')
	for l in lst:
		f.write("%s\n"%l)
	f.close()
	#np.savetxt('miriadextract_' + ano + "_sunpos.dat", sunpos[:len(t)], fmt='%8.5f')

###raw calibration################
if needrawcal and (not skip):
	for p, key in zip(range(len(wantpols)), wantpols.keys()):
		rawcalpar = np.fromfile(rawpaths[key], dtype="complex64").reshape(nfreq, nant)
		data[p] = omni.apply_calpar(data[p], rawcalpar, calibrators[p].totalVisibilityId)

####manually saving data################
if not skip:
	for p, calibrator in zip(range(len(wantpols)), calibrators):
		data[p].astype('complex64').tofile(calibrator.dataPath)
	del(data)
####calibrate################
print FILENAME + " MSG: starting calibration."
for p, calibrator in zip(range(len(wantpols)), calibrators):
	calibrator.nTime = len(timing)
	calibrator.nFrequency = nfreq
	calibrator.removeDegeneracy = removedegen
	calibrator.removeAdditive = removeadditive
	if calibrator.removeAdditive:
		calibrator.removeAdditivePeriod = removeadditiveperiod
	calibrator.keepData = keep_binary_data
	calibrator.keepCalpar = keep_binary_calpar
	calibrator.convergePercent = converge_percent
	calibrator.maxIteration = max_iter
	calibrator.stepSize = step_size
	print calibrator.nTime, calibrator.nFrequency
#	calibrator.readyForCpp()
	calibrator.loglincal(calibrator.dataPath, verbose=True)

