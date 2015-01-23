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
import cPickle as pickle
FILENAME = "first_cal.py"
PI = np.pi
TPI = 2 * np.pi
print "#Omnical Version %s#"%omni.__version__


######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

ap.scripting.add_standard_options(o, cal=True, pol=True)
#o.add_option('-d', '--datatag', action = 'store', default = 'PSA128', help = 'tag name of this data set')
#o.add_option('-i', '--infopath', action = 'store', default = 'DOESNTEXIST', help = 'Redundantinfo file to read.')
o.add_option('--max', action = 'store', type = 'int', default = 5, help = 'Max number of iterations when removing bad antennas.')
#o.add_option('--add', action = 'store_true', help = 'whether to enable crosstalk removal')
#o.add_option('--nadd', action = 'store', type = 'int', default = -1, help = 'time steps w to remove additive term with. for running average its 2w + 1 sliding window.')
#o.add_option('--datapath', action = 'store', default = None, help = 'uv file or binary file folder')
o.add_option('--healthbar', action = 'store', type = 'float', default = 2, help = 'Health threshold (0-100) over which an antenna is marked bad. 2 by default.')
o.add_option('--suppress', action = 'store', type = 'float', default = 1, help = 'Amplitude of the gains for the bad antennas. Larger means more suppressed.')
o.add_option('-f', '--freq_range', action = 'store', default = '0_0', help = 'Frequency bin number range to use for fitting amp and delay seperated by underscore. 0_0 by default and will process all frequencies.')
o.add_option('-o', '--outputpath', action = 'store', default = "DONT_WRITE", help = 'Output folder. No output by default.')
o.add_option('-t', '--info_tag', action = 'store', default = "DEFAULT", help = 'Name tag for output redundantinfo file.')
#o.add_option('-k', '--skip', action = 'store_true', help = 'whether to skip data importing from uv')
#o.add_option('-u', '--newuv', action = 'store_true', help = 'whether to create new uv files with calibration applied')
o.add_option('--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--ampdelay', action = 'store_true', help = 'whether to print out amplitude and delay for Calfile usage.')
#o.add_option('--smooth', action = 'store_true', help = 'whether to smooth the calibration results over frequency.')
o.add_option('--plot', action = 'store_true', help = 'Whether to make plots in the end.')
#o.add_option('--crude', action = 'store_true', help = 'whether to apply crude calibration')
o.add_option('-e', '--tol', action = 'store', type = 'float', default = 1e-2, help = 'tolerance of antenna location deviation when computing unique baselines.')
o.add_option('--ba', action = 'store', default = '', help = 'bad antenna number indices seperated by commas')
o.add_option('--bu', action = 'store', default = '', help = 'bad unique baseline indicated by ant pairs (seperated by .) seperated by commas: 1.2,3.4,10.11')



opts,args = o.parse_args(sys.argv[1:])
#skip = opts.skip
#create_new_uvs = opts.newuv
overwrite = opts.overwrite
make_plots = opts.plot
print_ampdelay = opts.ampdelay
#smooth = opts.smooth
#ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
#dataano = opts.datatag#ano for existing data and lst.dat
#sourcepath = opts.datapath
oppath = os.path.expanduser(opts.outputpath)
info_tag = opts.info_tag
uvfiles = args
healthbar = opts.healthbar
bad_ant_suppress = opts.suppress
max_try = opts.max

try:
	badAntenna = [int(i) for i in opts.ba.split(',')]
except:
	badAntenna = []
try:
	if opts.bu != '':
		badUBLpair = [[int(j) for j in i.split('.')] for i in opts.bu.split(',')]
	else:
		badUBLpair = []
except:
	badUBLpair = []
redundancy_tol = opts.tol

[fstart,fend] = [int(x) for x in opts.freq_range.split('_')]
for uvf in uvfiles:
	if not os.path.isdir(uvf):
		uvfiles.remove(uvf)
		print "WARNING: uv file path %s does not exist!"%uvf
if len(uvfiles) == 0:
	raise Exception("ERROR: No valid uv files detected in input. Exiting!")

wantpols = {}
for p in opts.pol.split(','): wantpols[p] = ap.miriad.str2pol[p]
#wantpols = {'xx':ap.miriad.str2pol['xx']}#, 'yy':-6}#todo:

print "Reading calfile %s..."%opts.cal,
sys.stdout.flush()
aa = ap.cal.get_aa(opts.cal, np.array([.15]))
print "Done. Antenna layout:"
print aa.ant_layout
sys.stdout.flush()


#infopaths = {}
#for pol in wantpols.keys():
	#infopaths[pol]= opts.infopath


removedegen = True #this is only for amplitude, because for phase there's a special degenaracy removal anyways
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
if fend == 0:
	fstart = int(nfreq) / 4
	fend = 3 * int(nfreq) / 4

nant = uv['nants']
sa = ephem.Observer()
sa.lon = uv['longitu']
sa.lat = uv['latitud']
sa.pressure = 0
startfreq = uv['sfreq']
dfreq = uv['sdf']
del(uv)
print "Done."
sys.stdout.flush()




###start reading miriads################
print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed"
sys.stdout.flush()
rawdata, t, timing, lst = omni.importuvs(uvfiles, np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), wantpols, timingTolerance=100)#, nTotalAntenna = len(aa))
print FILENAME + " MSG:",  len(t), "slices read."
sys.stdout.flush()

sun = ephem.Sun()
sunpos  = np.zeros((len(timing), 2))
southern_points = {'hyd':{'ra': '09:18:05.7', 'dec': '-12:05:44'},
'cen':{'ra': '13:25:27.6', 'dec': '-43:01:09'},
'cyg':{'ra': '19:59:28.3', 'dec': '40:44:02'},
'pic':{'ra': '05:19:49.7', 'dec': '-45:46:44'},
'vir':{'ra': '12:30:49.4', 'dec': '12:23:28'},
'for':{'ra': '03:22:41.7', 'dec': '-37:12:30'}}

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
new_bad_ant = ["Just to get while loop started"]
trials = 0
calibrators = {}
data = {}
while new_bad_ant != [] and trials < max_try:
	trials = trials + 1
	if trials > 1:
		print "##########################################################################"
		print FILENAME + " trial #%i: Recalculating redundant info removing new bad antennas..."%trials, new_bad_ant
		sys.stdout.flush()


	ant_bad_meter = {}
	ubl_bad_meter = {}
	crude_calpar = {}
	if trials > 1:
		#figure out old and new bad antennas
		badAntenna = list(np.sort(badAntenna + new_bad_ant))
		print 'Current bad Antennas (%i):'%len(badAntenna), badAntenna
		print 'Bad unique baselines (%i):'%len(badUBLpair), badUBLpair

	for p, pol in enumerate(wantpols.keys()):


		calibrators[pol] = omni.RedundantCalibrator_PAPER(aa)
		calibrators[pol].nTime = len(timing)
		calibrators[pol].nFrequency = nfreq
		calibrators[pol].removeDegeneracy = removedegen
		calibrators[pol].convergePercent = converge_percent
		calibrators[pol].maxIteration = max_iter
		calibrators[pol].stepSize = step_size


		timer = time.time()
		calibrators[pol].compute_redundantinfo(badAntenna = badAntenna, badUBLpair = badUBLpair, antennaLocationTolerance = redundancy_tol)
		print "Redundant info on %s computed in %f minutes."%(pol, (time.time() - timer)/60.)
		info = calibrators[pol].Info.get_info()

		###prepare rawCalpar for each calibrator and consider, if needed, raw calibration################
		if need_crude_cal:
			initant, solution_path, additional_solution_path, degen, _ = omni.find_solution_path(info, tol = calibrators[pol].antennaLocationTolerance, verbose = False)
			crude_calpar[pol] = np.array([omni.raw_calibrate(rawdata[p, 0, f], info, initant, solution_path, additional_solution_path, degen) for f in range(calibrators[pol].nFrequency)])
			data[p] = omni.apply_calpar(rawdata[p], crude_calpar[pol], calibrators[pol].totalVisibilityId)
		else:
			data[p] = rawdata[p]
		#calibrators[pol].rawCalpar = np.zeros((calibrators[pol].nTime, calibrators[pol].nFrequency, 3 + 2 * (calibrators[pol].Info.nAntenna + calibrators[pol].Info.nUBL)),dtype='float32')

		################################
		########calibrate################
		################################

		################first round of calibration	#########################
		print FILENAME + " MSG: starting calibration on %s. nTime = %i, nFrequency = %i ..."%(pol, calibrators[pol].nTime, calibrators[pol].nFrequency),
		sys.stdout.flush()
		timer = time.time()
		additivein = np.zeros_like(data[p])
		calibrators[pol].logcal(data[p], additivein, verbose=True)
		additiveout = calibrators[pol].lincal(data[p], additivein, verbose=True)
		print "Done. %fmin"%(float(time.time()-timer)/60.)
		sys.stdout.flush()

		#################degeneracy removal on 3 antennas [initant, degen[0], degen[1]]######################################
		print FILENAME + " MSG: Performing additional degeneracy removal on %s"%pol,
		sys.stdout.flush()
		timer = time.time()
		A = np.zeros((calibrators[pol].Info.nAntenna, 4))
		masker = np.zeros((calibrators[pol].Info.nAntenna, calibrators[pol].Info.nAntenna))
		for a in [initant, degen[0], degen[1]]:
			A[a] = list(calibrators[pol].Info.antloc[a]) + [1.]
			masker[a,a] = 1.
		matrix = np.identity(calibrators[pol].Info.nAntenna) - np.array([list(vec) + [1.] for vec in info['antloc']]).dot(la.pinv(A.transpose().dot(A)).dot(A.transpose()).dot(masker))

		calibrators[pol].rawCalpar[:, :, (3 + calibrators[pol].Info.nAntenna):(3 + 2 * calibrators[pol].Info.nAntenna)] = ((matrix.dot(calibrators[pol].rawCalpar[:, :, (3 + calibrators[pol].Info.nAntenna):(3 + 2 * calibrators[pol].Info.nAntenna)].transpose(0,2,1)).transpose(1,2,0)) + PI)%TPI - PI
		print "Done. %fmin"%(float(time.time()-timer)/60.)
		sys.stdout.flush()

		#######################diagnose###############################
		ant_bad_meter[pol], ubl_bad_meter[pol] = calibrators[pol].diagnose(data = data[p], additiveout = additiveout, healthbar = healthbar, verbose = False)
		nbad = 0
		for ab in ant_bad_meter[pol]:
			if ab > healthbar:
				nbad += 1
		if nbad > 0:
			print FILENAME + " MSG: %i bad antennas found on %s:"%(nbad, pol),
			for i, ab in enumerate(ant_bad_meter[pol]):
				if ab > healthbar:
					print calibrators[pol].Info.subsetant[i],
			print ""
		else:
			print FILENAME + " MSG: %i bad antennas found on %s"%(nbad, pol)
		sys.stdout.flush()


	new_bad_ant = []
	for a in range(calibrators[wantpols.keys()[0]].Info.nAntenna):
		for pol in wantpols.keys():
			if ant_bad_meter[pol][a] > healthbar:
				new_bad_ant.append(calibrators[wantpols.keys()[0]].Info.subsetant[a])
				break

print "Identified possible bad baselines (not automatically excluded):"
for u in range(calibrators[wantpols.keys()[0]].nUBL):#assuming calibrators on all pols are the same ubl config
	for p, pol in enumerate(wantpols.keys()):
		c = calibrators[pol]
		if ubl_bad_meter[pol][u] > 1:
			print np.round(c.ubl[u]), ubl_bad_meter[pol][u], "ant pair %i.%i"%(c.subsetant[int(c.ublindex[u][0][0])], c.subsetant[int(c.ublindex[u][0][1])])
			break

linearcalpar = {}#combine lincal results averaged over time and the initial crude_calpar results
for pol in wantpols.keys():
	linearcalpar[pol] = 10.**(calibrators[pol].rawCalpar[:,:,3:3+calibrators[pol].Info.nAntenna]) * np.exp(1j * calibrators[pol].rawCalpar[:,:,3+calibrators[pol].Info.nAntenna:3+2*calibrators[pol].Info.nAntenna])
	linearcalpar[pol] = nanmedian(np.real(linearcalpar[pol]), axis = 0) + 1j * nanmedian(np.imag(linearcalpar[pol]), axis = 0)
	linearcalpar[pol][np.isnan(linearcalpar[pol])] = 1
	linearcalpar[pol] = linearcalpar[pol] * crude_calpar[pol][:, calibrators[pol].Info.subsetant]

	#try smoothing out this calpar so it doesnt introduce much frequency structure

	freq_flag = (np.sum(calibrators[pol].flag(), axis = 0) > .4 * calibrators[pol].nTime)|np.isnan(np.sum(linearcalpar[pol], axis=1))
	A = np.array([list(a) + [1] for a in calibrators[pol].antloc])
	AAA = A.dot(np.linalg.pinv(A.transpose().dot(A)).dot(A.transpose()))
	#if smooth:
		#calpar = linearcalpar[pol][~freq_flag]
		#for i in range(15):
			#for f in range(len(calpar) - 1):
				#calpar[f+1] = np.exp(1.j*(np.angle(calpar[f+1]) - AAA.dot((np.angle(calpar[f+1]) - np.angle(calpar[f]) + PI)%TPI-PI)))
			#linearcalpar[pol][~freq_flag] = calpar
			#for f in np.array(range(calibrators[pol].nFrequency))[freq_flag]-1:
				#linearcalpar[pol][f+1] = np.exp(1.j*(np.angle(linearcalpar[pol][f+1]) - AAA.dot((np.angle(linearcalpar[pol][f+1]) - np.angle(linearcalpar[pol][f]) + PI)%TPI-PI)))


	####amplitude: a single amplitude over all frequency is only useful for calfile and not really accurate
	if print_ampdelay:
		amp = np.ones(calibrators[pol].nTotalAnt, dtype='float') * bad_ant_suppress
		amp[calibrators[pol].Info.subsetant] = nanmedian(np.abs(linearcalpar[pol]), axis = 0)
		print FILENAME + " MSG: amplitude factor on %s as |g|:"%pol
		print '{'
		for a1, a2 in zip(range(len(amp)), amp):
			print "%i: %f, "%(a1,a2)
		print '}'
		sys.stdout.flush()

	count = 0
	error_history = []
	delay = np.zeros(calibrators[pol].nTotalAnt, dtype='float')
	delay_error = np.zeros(calibrators[pol].nTotalAnt, dtype='float')+np.inf
	while (count == 0) or (np.max(delay_error[calibrators[pol].Info.subsetant]) > .1 and count < 100):
		count = count + 1
		####delay


		sub_freq_flag = freq_flag[fstart:fend]
		avg_angle = np.angle(linearcalpar[pol][fstart:fend].transpose()).astype('float32')#2D nant x freq

		####unwrap phase
		valid_freq_index = np.arange(len(freq_flag))[fstart:fend][~sub_freq_flag]
		rough_slopes = np.median(((avg_angle[:, ~sub_freq_flag][:, 1:] - avg_angle[:, ~sub_freq_flag][:, :-1] + PI)%TPI - PI) / (valid_freq_index[None, 1:]-valid_freq_index[None, :-1]), axis = 1)
		line_model = np.outer(rough_slopes, valid_freq_index)
		avg_angle[:, ~sub_freq_flag] = avg_angle[:, ~sub_freq_flag] + TPI * np.round((line_model - avg_angle[:, ~sub_freq_flag]) / TPI)
		#for a in range(len(avg_angle)):
			#avg_angle[a][~sub_freq_flag] = _O.unwrap_phase(avg_angle[a][~sub_freq_flag])


		##I am allowing a phase offset per antenna. iterating on degeneracy removal does not seem to help
		A = np.ones((fend-fstart, 2),dtype='float32')
		A[:, 0] = np.arange(fstart, fend) * dfreq + startfreq
		A = A[~sub_freq_flag]
		#####matrix_f0 = (la.pinv(A.transpose().dot(A)).dot(A.transpose()))[1]
		#####offsets = np.array([matrix_f0.dot(x[~sub_freq_flag]) for x in avg_angle])

		#####intersect = np.round(offsets/(2.*PI)) * (2.*PI) #find closest multiple of 2pi for the intersect
		#####overall_angle = 0#omni.meanAngle(offsets, weights = 1/(intersect+1))
		#####avg_angle = avg_angle - intersect[:,None] - overall_angle

		#######now fit
		#####A = np.arange(fstart, fend) * dfreq + startfreq
		#####A = A[~sub_freq_flag].reshape((np.sum(~sub_freq_flag), 1))
		matrix = (la.pinv(A.transpose().dot(A)).dot(A.transpose()))

		solution = matrix.dot(avg_angle[:, ~sub_freq_flag].transpose())
		delay[calibrators[pol].Info.subsetant] = solution[0] / TPI

		##use delay fit to further correct rephasing degeneracies in linearcalpar
		avg_angle_error = np.angle(linearcalpar[pol]) - np.outer(np.arange(nfreq) * dfreq + startfreq, solution[0]) - solution[1]
		avg_angle_error = (avg_angle_error + PI)%TPI - PI
		fit_degen = avg_angle_error.dot(AAA.transpose())
		fit_degen[np.isnan(fit_degen)] = 0
		linearcalpar[pol] = linearcalpar[pol] / np.exp(1.j*fit_degen)
		avg_angle_error = (avg_angle_error - fit_degen + PI)%TPI - PI
		delay_error[calibrators[pol].Info.subsetant] = np.linalg.norm(avg_angle_error[fstart:fend][~sub_freq_flag,:], axis = 0) / (np.sum(~sub_freq_flag))**.5
		error_history + []

	print count



	#print stuff for calfile
	if print_ampdelay:
		print FILENAME + " MSG: delay on %s in nanoseconds:"%pol
		print '{'
		for a1, a2 in zip(range(len(delay)), delay):
			print "%i: %f, "%(a1,a2)
		print '}'
		sys.stdout.flush()

	if make_plots:
		nplot = 8
		plot_a = np.argsort(delay_error[calibrators[pol].subsetant])[range(0, calibrators[pol].nAntenna - 1, calibrators[pol].nAntenna/min(nplot,calibrators[pol].nAntenna)) + [-1]]#[range(nplot/2)+range(-nplot/2,0)]
		for i, a in enumerate(plot_a):
			plt.subplot(1, len(plot_a), i+1)
			plt.plot(np.arange(nfreq), (np.angle(linearcalpar[pol][:, a]) + PI)%TPI - PI)
			plt.plot(np.arange(nfreq), ((np.arange(nfreq)*dfreq + startfreq) * solution[0, a] + solution[1, a] + PI)%TPI - PI)
			plt.title("Ant#%i %.1f"%(calibrators[pol].subsetant[a], delay_error[calibrators[pol].Info.subsetant[a]]))
			plt.axis([0, nfreq, -PI, PI])
			#plt.axes().set_aspect('equal')
		plt.show()
		plt.hist(delay_error[calibrators[pol].Info.subsetant], 20)
		plt.show()


if oppath != "DONT_WRITE/":
	if info_tag == "DEFAULT":
		op_info_path = oppath + 'redundantinfo_first_cal_' + time.strftime("%Y_%m_%d_%H_%M_%S") + ".bin"
		op_calpar_path = oppath + 'calpar_first_cal_' + time.strftime("%Y_%m_%d_%H_%M_%S") + ".p"
	else:
		op_info_path = oppath + 'redundantinfo_first_cal_' + info_tag + ".bin"
		op_calpar_path = oppath + 'calpar_first_cal_' + info_tag + ".p"

	print FILENAME + " MSG: Writing redundant info to %s"%op_info_path,
	sys.stdout.flush()
	calibrators[wantpols.keys()[0]].write_redundantinfo(infoPath = op_info_path, verbose = False, overwrite = overwrite)
	print "Done."
	sys.stdout.flush()

	#for pol in wantpols.keys():
		#print FILENAME + " MSG: Writing %s raw calpar to %s"%(pol, op_calpar_path.replace('.bin', pol+'.bin')),
		#sys.stdout.flush()
		#linearcalpar[pol].astype('float32').tofile(op_calpar_path.replace('.bin', pol+'.bin'))
		#print "Done."
		#sys.stdout.flush()
	linearcalpar_out = {}#include all antennas, not just good ones
	for pol in linearcalpar.keys():
		linearcalpar_out[pol] = np.ones((linearcalpar[pol].shape[0], calibrators[pol].nTotalAnt), dtype='complex64')
		linearcalpar_out[pol][:, calibrators[pol].subsetant] = linearcalpar[pol]
		linearcalpar_out[pol][np.isnan(linearcalpar_out[pol])] = 1
	with open(op_calpar_path, 'wb') as outfile:

		pickle.dump(linearcalpar_out, outfile, protocol=pickle.HIGHEST_PROTOCOL)
else:
	print FILENAME + " MSG: Not outputting redundantinfo or rawcalpar by default."
	sys.stdout.flush()





print "%i Bad antennas found: "%len(badAntenna),
bad_str = ""
for bant in badAntenna:
	bad_str += (str(int(bant)) + ",")
if bad_str[-1] == ",":
	bad_str = bad_str[:-1]
print bad_str
