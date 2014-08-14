#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys
import scipy.signal as ss
import scipy.linalg as la
from scipy.stats import nanmedian
import matplotlib.pyplot as plt
FILENAME = "first_cal.py"

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
#o.add_option('-t', '--tag', action = 'store', default = 'PSA128', help = 'tag name of this calibration')
#o.add_option('-d', '--datatag', action = 'store', default = 'PSA128', help = 'tag name of this data set')
o.add_option('-i', '--infopath', action = 'store', default = '/data2/home/hz2ug/omnical/doc/redundantinfo_PSA128_17ba.bin', help = 'redundantinfo file to read')
#o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
#o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
#o.add_option('--datapath', action = 'store', default = None, help = 'uv file or binary file folder')
o.add_option('--healthbar', action = 'store', type = 'float', default = 10, help = 'health threshold (0-100) over which an antenna is marked bad.')
o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
#o.add_option('-k', '--skip', action = 'store_true', help = 'whether to skip data importing from uv')
#o.add_option('-u', '--newuv', action = 'store_true', help = 'whether to create new uv files with calibration applied')
#o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--plot', action = 'store_true', help = 'whether to make plots in the end')
#o.add_option('--crude', action = 'store_true', help = 'whether to apply crude calibration')

opts,args = o.parse_args(sys.argv[1:])
#skip = opts.skip
#create_new_uvs = opts.newuv
#overwrite_uvs = opts.overwrite
make_plots = opts.plot
#ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
#dataano = opts.datatag#ano for existing data and lst.dat
#sourcepath = opts.datapath
oppath = opts.outputpath
uvfiles = args
healthbar = opts.healthbar
for uvf in uvfiles:
	if not os.path.isdir(uvf):
		uvfiles.remove(uvf)
		print "WARNING: uv file path %s does not exist!"%uvf
if len(uvfiles) == 0:
	raise Exception("ERROR: No valid uv files detected in input. Exiting!")

wantpols = {}
for p in opts.pol.split(','): wantpols[p] = ap.miriad.str2pol[p]
#wantpols = {'xx':ap.miriad.str2pol['xx']}#, 'yy':-6}#todo:

print "Reading calfile %s"%opts.cal,
sys.stdout.flush()
aa = ap.cal.get_aa(opts.cal, np.array([.15]))
print "Done"
sys.stdout.flush()


infopaths = {}
for key in wantpols.keys():
	infopaths[key]= opts.infopath


removedegen = False
removeadditive = False
removeadditiveperiod = -1

need_crude_cal = True

converge_percent = 0.01
max_iter = 50
step_size = .3

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
oppath += '/'

####get some info from the first uvfile   ################
print "Getting some basic info from %s"%uvfiles[0],
sys.stdout.flush()
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants']
sa = ephem.Observer()
sa.lon = uv['longitu']
sa.lat = uv['latitud']
#startfreq = uv['sfreq']
dfreq = uv['sdf']
del(uv)
print "Done."
sys.stdout.flush()




###start reading miriads################
print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed"
sys.stdout.flush()
data, t, timing, lst = omni.importuvs(uvfiles, np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), wantpols, timingTolerance=100)#, nTotalAntenna = len(aa))
print FILENAME + " MSG:",  len(t), "slices read."
sys.stdout.flush()

sun = ephem.Sun()
sunpos  = np.zeros((len(timing), 2))
southern_points = {'hyd':{'ra': '09:18:05.7', 'dec': '-12:05:44'},
'cen':{'ra': '13:25:27.6', 'dec': '-43:01:09'},
'cyg':{'ra': '19:59:28.3', 'dec': '40:44:02'},
'pic':{'ra': '05:19:49.7', 'dec': '-45:46:44'},
'vir':{'ra': '12:30:49.4', 'dec': '12:23:28'},
'for':{'ra': '03:22:41.7', 'dec': '-37:12m30s'}}

for source in southern_points.keys():
	southern_points[source]['body'] = ephem.FixedBody()
	southern_points[source]['body']._ra = southern_points[source]['ra']
	southern_points[source]['body']._dec = southern_points[source]['dec']
	southern_points[source]['pos'] = np.zeros((len(timing), 2))
for nt,tm in zip(range(len(timing)),timing):
	sa.date = tm
	sun.compute(sa)
	sunpos[nt] = sun.alt, sun.az
	for source in southern_points.keys():
		southern_points[source]['body'].compute(sa)
		southern_points[source]['pos'][nt] = southern_points[source]['body'].alt, southern_points[source]['body'].az
print FILENAME + " MSG:"
print "data time range UTC: %s to %s"%(timing[0], timing[-1])
print "sun altaz from (%f,%f) to (%f,%f)"%(sunpos[0,0], sunpos[0,1], sunpos[-1,0], sunpos[-1,1])
for source in southern_points.keys():
	print "%s altaz from (%f,%f) to (%f,%f)"%(source, southern_points[source]['pos'][0,0], southern_points[source]['pos'][0,1], southern_points[source]['pos'][-1,0], southern_points[source]['pos'][-1,1])
sys.stdout.flush()
####create redundant calibrators################
#calibrators = [omni.RedundantCalibrator(nant, info = infopaths[key]) for key in wantpols.keys()]
calibrators = {}
ant_bad_meter = {}
crude_calpar = {}
for p, key in zip(range(len(data)), wantpols.keys()):

	calibrators[key] = RedundantCalibrator_PAPER(aa)
	calibrators[key].read_redundantinfo(infopaths[key], verbose=False)
	info = calibrators[key].Info.get_info()
	calibrators[key].nTime = len(timing)
	calibrators[key].nFrequency = nfreq

	###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
	if need_crude_cal:
		initant, solution_path, additional_solution_path, degen, _ = omni.find_solution_path(info)
		crude_calpar[key] = np.array([omni.raw_calibrate(data[p, 0, f], info, initant, solution_path, additional_solution_path, degen) for f in range(calibrators[key].nFrequency)])
		data[p] = omni.apply_calpar(data[p], crude_calpar[key], calibrators[key].totalVisibilityId)

	calibrators[key].rawCalpar = np.zeros((calibrators[key].nTime, calibrators[key].nFrequency, 3 + 2 * (calibrators[key].Info.nAntenna + calibrators[key].Info.nUBL)),dtype='float32')
	####calibrate################
	calibrators[key].removeDegeneracy = removedegen
	calibrators[key].convergePercent = converge_percent
	calibrators[key].maxIteration = max_iter
	calibrators[key].stepSize = step_size

	################first round of calibration	#########################
	print FILENAME + " MSG: starting calibration on %s. nTime = %i, nFrequency = %i ..."%(key, calibrators[key].nTime, calibrators[key].nFrequency),
	sys.stdout.flush()
	timer = time.time()
	additivein = np.zeros_like(data[p])
	calibrators[key].logcal(data[p], additivein, verbose=True)
	additiveout = calibrators[key].lincal(data[p], additivein, verbose=True)
	print "Done. %fmin"%(float(time.time()-timer)/60.)
	sys.stdout.flush()
	#######################diagnose###############################
	ant_bad_meter[key], _ = calibrators[key].diagnose(data = data[p], additiveout = additiveout, healthbar = healthbar, verbose = False)
	nbad = 0
	for ab in ant_bad_meter[key]:
		if ab > healthbar:
			nbad += 1
	print FILENAME + " MSG: %i badantennas found on %s"%(nbad, key)
	sys.stdout.flush()


new_bad_ant = []
for a in range(calibrators[wantpols.keys()[0]].Info.nAntenna):
	for key in wantpols.keys():
		if ant_bad_meter[key][a] > healthbar:
			new_bad_ant.append(calibrators[wantpols.keys()[0]].Info.subsetant[a])
			break
if new_bad_ant != []:
	print FILENAME + " MSG: Recalculating redundant info removing badantennas..."
	sys.stdout.flush()
	print "Done."
	sys.stdout.flush()

####amplitude
for p,pol in zip(range(len(wantpols)), wantpols.keys()):
	amp = np.ones(calibrators[pol].nTotalAnt, dtype='float')
	amp[calibrators[pol].Info.subsetant] = 10**(nanmedian(nanmedian(calibrators[pol].rawCalpar[:,:,3:3+calibrators[pol].Info.nAntenna],axis=0),axis=0))
	print FILENAME + " MSG: amplitude factor on %s:"%pol, amp
	sys.stdout.flush()

####delay
for p,pol in zip(range(len(wantpols)), wantpols.keys()):
	delay = np.zeros(calibrators[pol].nTotalAnt, dtype='float')
	delay_error = np.zeros(calibrators[pol].nTotalAnt, dtype='float')+np.inf
	fstart = 80
	fend = 150
	A = np.ones((fend-fstart, 2),dtype='float32')
	A[:, 0] = range(fstart, fend)
	matrix = (la.pinv(A.transpose().dot(A)).dot(A.transpose()))[0]
	error_matrix = A.dot(la.pinv(A.transpose().dot(A)).dot(A.transpose())) - np.identity(len(A))
	avg_angle = np.angle((np.nanmean(np.exp(1.j * calibrators[pol].rawCalpar[:,:,3+calibrators[pol].Info.nAntenna:3+2*calibrators[pol].Info.nAntenna]), axis = 0) * crude_calpar[key][:, calibrators[pol].Info.subsetant])[fstart:fend].transpose())#2D nant x freq
	avg_angle -= avg_angle[0]
	nplot = 10
	if make_plots:
		for a in range(0, len(avg_angle), len(avg_angle)/min(10,len(avg_angle))):
			plt.subplot(1, min(10,len(avg_angle)), (a/( len(avg_angle)/min(10,len(avg_angle)))))
			plt.plot(avg_angle[a])
			plt.plot((np.arange(len(A)) * matrix.dot(_O.unwrap_phase( avg_angle[a])) + np.pi)%(2*np.pi) - np.pi)
			plt.axis([0, len(A), -np.pi, np.pi])
		plt.show()
	delay[calibrators[pol].Info.subsetant] = [matrix.dot(_O.unwrap_phase(x))/ (2 * np.pi * dfreq)  for x in avg_angle]
	delay_error[calibrators[pol].Info.subsetant] = [la.norm(error_matrix.dot(_O.unwrap_phase(x)))/ (len(A))**.5 for x in avg_angle]
	print FILENAME + " MSG: delay on %s:"%pol
	print delay
	print delay_error

	sys.stdout.flush()
#if make_plots:
	#for p,pol in zip(range(len(wantpols)), wantpols.keys()):
		#plt.subplot(1,len(wantpols),p+1)
		#plot_data = (calibrators[pol].rawCalpar[:,:,2]/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5
		#plt.imshow(plot_data, vmin = 0, vmax = (np.nanmax(calibrators[wantpols.keys()[0]].rawCalpar[:,30:-30:5,2])/(len(calibrators[pol].Info.subsetbl)-calibrators[pol].Info.nAntenna - calibrators[pol].Info.nUBL))**.5, interpolation='nearest')
	#plt.title('RMS fitting error per baseline')
	#plt.colorbar()
	#plt.show()
