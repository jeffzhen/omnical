import datetime
import socket, multiprocessing, math, random, traceback, ephem, string, commands, datetime, shutil
import time
from time import ctime
import aipy as ap
import struct
import numpy as np
import os, sys
from optparse import OptionParser
import omnical._omnical as _O
import warnings
from array import array
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy as sp
    import scipy.sparse as sps
    import scipy.linalg as la

FILENAME = "calibration_omni.py"


infokeys = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']

int_infokeys = ['nAntenna','nUBL','nBaseline']
intarray_infokeys = ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
float_infokeys = ['antloc','ubl','degenM','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
def read_redundantinfo(infopath, verbose = False):
	METHODNAME = "read_redundantinfo"
	if not os.path.isfile(infopath):
		raise Exception('Error: file path %s does not exist!'%infopath)
	timer = time.time()
	with open(infopath) as f:
		rawinfo = np.array([np.array([float(x) for x in line.split()]) for line in f])
	if len(rawinfo) < len(infokeys):
		raise Exception('Error: number of rows in %s (%i) is less than expected length of %i!'%(infopath, len(rawinfo), len(infokeys)))
	if verbose:
		print FILENAME + "*" + METHODNAME + " MSG:",  "Reading redundant info...",

	info = {}
	infocount = 0;
	info['nAntenna'] = int(rawinfo[infocount][0]) #number of good antennas among all (64) antennas, same as the length of subsetant
	infocount += 1

	info['nUBL'] = int(rawinfo[infocount][0]) #number of unique baselines
	infocount += 1

	nbl = int(rawinfo[infocount][0])
	info['nBaseline'] = nbl
	infocount += 1


	info['subsetant'] = rawinfo[infocount].astype(int) #the index of good antennas in all (64) antennas
	infocount += 1

	info['antloc'] = rawinfo[infocount].reshape((info['nAntenna'],3)) #the index of good antennas in all (64) antennas
	infocount += 1

	info['subsetbl'] = rawinfo[infocount].astype(int) #the index of good baselines (auto included) in all baselines
	infocount += 1
	info['ubl'] = rawinfo[infocount].reshape((info['nUBL'],3)) #unique baseline vectors
	infocount += 1
	info['bltoubl'] = rawinfo[infocount].astype(int) #cross bl number to ubl index
	infocount += 1
	info['reversed'] = rawinfo[infocount].astype(int) #cross only bl if reversed -1, otherwise 1
	infocount += 1
	info['reversedauto'] = rawinfo[infocount].astype(int) #the index of good baselines (auto included) in all baselines
	infocount += 1
	info['autoindex'] = rawinfo[infocount].astype(int)  #index of auto bls among good bls
	infocount += 1
	info['crossindex'] = rawinfo[infocount].astype(int)  #index of cross bls among good bls
	infocount += 1
	ncross = len(info['crossindex'])
	#info['ncross'] = ncross
	info['bl2d'] = rawinfo[infocount].reshape(nbl, 2).astype(int) #from 1d bl index to a pair of antenna numbers
	infocount += 1
	info['ublcount'] = rawinfo[infocount].astype(int) #for each ubl, the number of good cross bls corresponding to it
	infocount += 1
	info['ublindex'] = range((info['nUBL'])) #//for each ubl, the vector<int> contains (ant1, ant2, crossbl)
	tmp = rawinfo[infocount].reshape(ncross, 3).astype(int)
	infocount += 1
	cnter = 0
	for i in range(info['nUBL']):
		info['ublindex'][i] = np.zeros((info['ublcount'][i],3))
		for j in range(len(info['ublindex'][i])):
			info['ublindex'][i][j] = tmp[cnter]
			cnter+=1
	info['ublindex'] = np.asarray(info['ublindex'])

	info['bl1dmatrix'] = rawinfo[infocount].reshape((info['nAntenna'], info['nAntenna'])).astype(int) #a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
	infocount += 1
	#matrices
	info['degenM'] = rawinfo[infocount].reshape((info['nAntenna'] + info['nUBL'], info['nAntenna']))
	infocount += 1
	info['A'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #A matrix for logcal amplitude
	infocount += 1
	info['B'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #B matrix for logcal phase
	infocount += 1
	##The sparse matrices are treated a little differently because they are not rectangular
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore",category=DeprecationWarning)
		info['At'] = info['A'].transpose()
		info['Bt'] = info['B'].transpose()
		info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense(), cond = 10**(-6))#(AtA)^-1
		info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense(), cond = 10**(-6))#(BtB)^-1
		info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
		info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
		info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
		info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
		info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
		info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB
	if verbose:
		print "done. nAntenna, nUBL, nBaseline = %i, %i, %i. Time taken: %f minutes."%(len(info['subsetant']), info['nUBL'], info['nBaseline'], (time.time()-timer)/60.)
	return info


def write_redundantinfo(info, infopath, overwrite = False, verbose = False):
	METHODNAME = "*write_redundantinfo*"
	timer = time.time()
	if (not overwrite) and os.path.isfile(infopath):
		raise Exception("Error: a file exists at " + infopath + ". Use overwrite = True to overwrite.")
		return
	if (overwrite) and os.path.isfile(infopath):
		os.remove(infopath)
	f_handle = open(infopath,'a')
	for key in infokeys:
		if key in ['antloc', 'ubl', 'degenM', 'AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:
			np.savetxt(f_handle, [np.array(info[key]).flatten()])
		elif key == 'ublindex':
			np.savetxt(f_handle, [np.vstack(info[key]).flatten()], fmt = '%d')
		elif key in ['At','Bt']:
			tmp = []
			for i in range(info[key].shape[0]):
				for j in range(info[key].shape[1]):
					if info[key][i,j] != 0:
						tmp += [i, j, info[key][i,j]]
			np.savetxt(f_handle, [np.array(tmp).flatten()], fmt = '%d')
		elif key in ['A','B']:
			np.savetxt(f_handle, info[key].todense().flatten(), fmt = '%d')
		else:
			np.savetxt(f_handle, [np.array(info[key]).flatten()], fmt = '%d')
	f_handle.close()
	if verbose:
		print "Info file successfully written to %s. Time taken: %f minutes."%(infopath, (time.time()-timer)/60.)
	return

def write_redundantinfo_binary(info, infopath, overwrite = False, verbose = False):
	METHODNAME = "*write_redundantinfo*"
	timer = time.time()
	if (not overwrite) and os.path.isfile(infopath):
		raise Exception("Error: a file exists at " + infopath + ". Use overwrite = True to overwrite.")
		return
	if (overwrite) and os.path.isfile(infopath):
		os.remove(infopath)
	outfile=open(infopath,'wb')
	marker = 9999999
	array('d',[marker]).tofile(outfile)     #start with a marker
	for key in infokeys:
		if key in ['antloc', 'ubl','degenM', 'AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:  #'antloc',
			farray = array('d',np.array(info[key]).flatten())
			farray.tofile(outfile)
			array('d',[marker]).tofile(outfile)
		elif key == 'ublindex':
			farray = array('d',np.vstack(info[key]).flatten())
			farray.tofile(outfile)
			array('d',[marker]).tofile(outfile)
		elif key in ['At','Bt']:
			tmp = []
			for i in range(info[key].shape[0]):
				for j in range(info[key].shape[1]):
					if info[key][i,j] != 0:
						tmp += [i, j, info[key][i,j]]
			farray = array('d',np.array(tmp).flatten())
			farray.tofile(outfile)
			array('d',[marker]).tofile(outfile)
		elif key in ['A','B']:
			farray = array('d',np.array(info[key].todense().flatten()).flatten())
			farray.tofile(outfile)
			array('d',[marker]).tofile(outfile)
		else:
			farray = array('d',np.array(info[key]).flatten())
			farray.tofile(outfile)
			array('d',[marker]).tofile(outfile)	
	outfile.close()
	if verbose:
		print "Info file successfully written to %s. Time taken: %f minutes."%(infopath, (time.time()-timer)/60.)
	return
	
def read_redundantinfo_binary(infopath, verbose = False):
	timer = time.time()
	METHODNAME = "read_redundantinfo"
	
	if not os.path.isfile(infopath):
		raise Exception('Error: file path %s does not exist!'%infopath)
	if verbose:
		print FILENAME + "*" + METHODNAME + " MSG:",  "Reading redundant info...",
	with open(infopath) as f:
		farray=array('d')
		farray.fromstring(f.read())
		datachunk = np.array(farray)
		marker = 9999999
		markerindex=np.where(datachunk == marker)[0]
		rawinfo=np.array([np.array(datachunk[markerindex[i]+1:markerindex[i+1]]) for i in range(len(markerindex)-1)])


	info = {}
	infocount = 0;
	info['nAntenna'] = int(rawinfo[infocount][0]) #number of good antennas among all (64) antennas, same as the length of subsetant
	infocount += 1
	info['nUBL'] = int(rawinfo[infocount][0]) #number of unique baselines
	infocount += 1
	nbl = int(rawinfo[infocount][0])
	info['nBaseline'] = nbl
	infocount += 1
	info['subsetant'] = rawinfo[infocount].astype(int) #the index of good antennas in all (64) antennas
	infocount += 1
	info['antloc'] = rawinfo[infocount].reshape((info['nAntenna'],3)) #the index of good antennas in all (64) antennas
	infocount += 1
	info['subsetbl'] = rawinfo[infocount].astype(int) #the index of good baselines (auto included) in all baselines
	infocount += 1
	info['ubl'] = rawinfo[infocount].reshape((info['nUBL'],3)) #unique baseline vectors
	infocount += 1
	info['bltoubl'] = rawinfo[infocount].astype(int) #cross bl number to ubl index
	infocount += 1
	info['reversed'] = rawinfo[infocount].astype(int) #cross only bl if reversed -1, otherwise 1
	infocount += 1
	info['reversedauto'] = rawinfo[infocount].astype(int) #the index of good baselines (auto included) in all baselines
	infocount += 1
	info['autoindex'] = rawinfo[infocount].astype(int)  #index of auto bls among good bls
	infocount += 1
	info['crossindex'] = rawinfo[infocount].astype(int)  #index of cross bls among good bls
	infocount += 1
	ncross = len(info['crossindex'])
	info['ncross'] = ncross
	info['bl2d'] = rawinfo[infocount].reshape(nbl, 2).astype(int) #from 1d bl index to a pair of antenna numbers
	infocount += 1
	info['ublcount'] = rawinfo[infocount].astype(int) #for each ubl, the number of good cross bls corresponding to it
	infocount += 1
	info['ublindex'] = range((info['nUBL'])) #//for each ubl, the vector<int> contains (ant1, ant2, crossbl)
	tmp = rawinfo[infocount].reshape(ncross, 3).astype(int)
	infocount += 1
	cnter = 0
	for i in range(info['nUBL']):
		info['ublindex'][i] = np.zeros((info['ublcount'][i],3))
		for j in range(len(info['ublindex'][i])):
			info['ublindex'][i][j] = tmp[cnter]
			cnter+=1
	info['ublindex'] = np.asarray(info['ublindex'])

	info['bl1dmatrix'] = rawinfo[infocount].reshape((info['nAntenna'], info['nAntenna'])).astype(int) #a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
	infocount += 1
	#matrices
	info['degenM'] = rawinfo[infocount].reshape((info['nAntenna'] + info['nUBL'], info['nAntenna']))
	infocount += 1
	info['A'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #A matrix for logcal amplitude
	infocount += 1
	info['B'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #B matrix for logcal phase
	infocount += 1
	##The sparse matrices are treated a little differently because they are not rectangular
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore",category=DeprecationWarning)
		info['At'] = info['A'].transpose()
		info['Bt'] = info['B'].transpose()
		info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense(), cond = 10**(-6))#(AtA)^-1
		info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense(), cond = 10**(-6))#(BtB)^-1
		info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
		info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
		info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
		info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
		info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
		info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB
	if verbose:
		print "done. nAntenna, nUBL, nBaseline = %i, %i, %i. Time taken: %f minutes."%(len(info['subsetant']), info['nUBL'], info['nBaseline'], (time.time()-timer)/60.)
	return info

def importuvs(uvfilenames, totalVisibilityId, wantpols, nTotalAntenna = None, timingTolerance = math.pi/12/3600/100):#tolerance of timing in radians in lst
	METHODNAME = "*importuvs*"
	############################################################
	sun = ephem.Sun()
	julDelta = 2415020 # =julian date - pyephem's Observer date
	####get some info from the first uvfile####################
	uv=ap.miriad.UV(uvfilenames[0])
	nfreq = uv.nchan;
	if nTotalAntenna == None:
		nant = uv['nants'] # 'nants' should be the number of dual-pol antennas. PSA32 has a bug in double counting
	else:
		nant = nTotalAntenna
	if nant * (nant + 1) / 2 < len(totalVisibilityId):
		raise Exception("FATAL ERROR: Total number of antenna %d implies %d baselines whereas the length of totalVisibilityId is %d."%(nant, nant * (nant + 1) / 2, len(totalVisibilityId)))
	startfreq = uv['sfreq']
	dfreq = uv['sdf']

	sa = ephem.Observer()
	sa.lon = uv['longitu']
	sa.lat = uv['latitud']

	del(uv)

	#######compute bl1dmatrix####each entry is 1 indexed with minus meaning conjugate
	bl1dmatrix = np.zeros((nant, nant), dtype = 'int32')
	for a1a2, bl in zip(totalVisibilityId, range(len(totalVisibilityId))):
		a1, a2 = a1a2
		bl1dmatrix[a1, a2] = bl + 1
		bl1dmatrix[a2, a1] = - (bl + 1)
	####prepare processing
	deftime = 2000
	data = np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
	#sunpos = np.zeros((deftime, 2))
	t = []#julian date
	timing = []#local time string
	lst = []#in units of sidereal hour

	###start processing
	datapulled = False
	for uvfile in uvfilenames:
		uv = ap.miriad.UV(uvfile)
		if len(timing) > 0:
			print FILENAME + METHODNAME + "MSG:",  timing[-1]#uv.nchan
		#print FILENAME + " MSG:",  uv['nants']
		currentpol = 0
		for preamble, rawd in uv.all():
			if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
				t += [preamble[1]]
				sa.date = preamble[1] - julDelta
				#sun.compute(sa)
				timing += [sa.date.__str__()]
				if abs((uv['lst'] - float(sa.sidereal_time()) + math.pi)%(2*math.pi) - math.pi) >= timingTolerance:
					raise Exception("Error: uv['lst'] is %f radians whereas time computed by aipy is %f radians, the difference is larger than tolerance of %f."%(uv['lst'], float(sa.sidereal_time()), timingTolerance))
				else:
					lst += [(float(sa.sidereal_time()) * 24./2./math.pi)]
				if len(t) > len(data):
					print FILENAME + METHODNAME + " MSG:",  "expanding number of time slices from", len(data), "to", len(data) + deftime
					data = np.concatenate((data, np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64')))
					#sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
					#sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
			for p, pol in zip(range(len(wantpols)), wantpols.keys()):
				if wantpols[pol] == uv['pol']:#//todo: use select()
					a1, a2 = preamble[2]
					bl = bl1dmatrix[a1, a2]#bl is 1 indexed with minus meaning conjugate
					datapulled = True
					#print info[p]['subsetbl'][info[p]['crossindex'][bl]],
					data[len(t) - 1, p, abs(bl) - 1] = (np.real(rawd.data) + 1.j * np.sign(bl) * np.imag(rawd.data)).astype('complex64')
		del(uv)
		if not datapulled:
			print FILENAME + METHODNAME + " MSG:",  "FATAL ERROR: no data pulled from " + uvfile + ", check polarization information! Exiting."
			exit(1)
	reorder = (1, 0, 3, 2)
	return np.transpose(data[:len(t)],reorder), t, timing, lst

def apply_calpar(data, calpar, visibilityID):#apply complex calpar for all antennas onto all baselines, calpar's dimension will be assumed to mean: 1D: constant over time and freq; 2D: constant over time; 3D: change over time and freq
	METHODNAME = "*apply_calpar*"
	if calpar.shape[-1] <= np.amax(visibilityID) or data.shape[-1] != len(visibilityID):
		raise Exception("Dimension mismatch! Either number of antennas in calpar " + str(calpar.shape[-1]) + " is less than implied in visibility ID "  + str(1 + np.amax(visibilityID)) + ", or the length of the last axis of data "  + str(data.shape[-1]) + " is not equal to length of visibilityID "  + str(len(visibilityID)) + ".")
	if len(calpar.shape) == 3 and len(data.shape) == 3 and calpar.shape[:2] == data.shape[:2]:
		return data/(np.conjugate(calpar[:,:,visibilityID[:,0].astype(int)]) * calpar[:,:,visibilityID[:,1].astype(int)])
	elif len(calpar.shape) == 2 and (len(data.shape) == 3 or len(data.shape) == 2) and calpar.shape[0] == data.shape[-2]:
		return data/(np.conjugate(calpar[:,visibilityID[:,0].astype(int)]) * calpar[:,visibilityID[:,1].astype(int)])
	elif len(calpar.shape) == 1 and len(data.shape) <= 3:
		return data/(np.conjugate(calpar[visibilityID[:,0].astype(int)]) * calpar[visibilityID[:,1].astype(int)])
	else:
		raise Exception("Dimension mismatch! I don't know how to interpret data dimension of " + str(data.shape) + " and calpar dimension of " + str(calpar.shape) + ".")

def apply_omnical_uvs(uvfilenames, calparfilenames, totalVisibilityId, info, wantpols, oppath, ano, additivefilenames = None, nTotalAntenna = None, overwrite= False):
	METHODNAME = "*apply_omnical_uvs*"
	if len(additivefilenames) != len(calparfilenames) and additivefilenames != None:
		raise Exception("Error: additivefilenames and calparfilenames have different lengths!")
	if len(info) != len(calparfilenames):
		raise Exception("Error: info and calparfilenames have different lengths!")
	if additivefilenames == None:
		additivefilenames = ["iDontThinkYouHaveAFileCalledThis" for _ in calparfilenames]

	####get some info from the first uvfile
	uv=ap.miriad.UV(uvfilenames[0])
	nfreq = uv.nchan;
	if nTotalAntenna == None:
		nant = uv['nants'] # 'nants' should be the number of dual-pol antennas. PSA32 has a bug in double counting
	else:
		nant = nTotalAntenna

	if nant * (nant + 1) / 2 < len(totalVisibilityId):
		raise Exception("FATAL ERROR: Total number of antenna %d implies %d baselines whereas the length of totalVisibilityId is %d."%(nant, nant * (nant + 1) / 2, len(totalVisibilityId)))
	startfreq = uv['sfreq']
	dfreq = uv['sdf']
	del(uv)

	#######compute bl1dmatrix####each entry is 1 indexed with minus meaning conjugate, the bl here is not number in totalVisibilityId, but in info['subsetbl'], so it's different from bl1dmatrix in import_uvs method. it also has 2 pols
	bl1dmatrix = [np.zeros((nant, nant), dtype = 'int32') for p in range(len(info))]

	for a1a2, bl in zip(totalVisibilityId, range(len(totalVisibilityId))):
		a1, a2 = a1a2
		for p in range(len(info)):
			for sbl, bl2 in zip(range(len(info[p]['subsetbl'])), info[p]['subsetbl']):
				if bl == bl2:
					bl1dmatrix[p][a1, a2] = sbl + 1
					bl1dmatrix[p][a2, a1] = - (sbl + 1)
					break
	####load calpar and check dimensions, massage calpar from txfx(3+2a+2u) to t*goodabl*f
	calpars = []#bad antenna included
	adds = []#badubl not included
	for p in range(len(wantpols)):
		calpar = np.fromfile(calparfilenames[p], dtype='float32')
		if len(calpar)%(nfreq *( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']))) != 0:
			print FILENAME + METHODNAME + " MSG:",  "FATAL ERROR: calpar input array " + calparfilenames[p] + " has length", calpar.shape, "which is not divisible by ", nfreq, 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']), "Aborted!"
			return
		ttotal = len(calpar)/(nfreq *( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL'])))
		calpar = calpar.reshape((ttotal, nfreq, ( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']))))
		calpars.append(1 + np.zeros((ttotal, nant, nfreq),dtype='complex64'))
		calpars[p][:,info[p]['subsetant'],:] = ((10**calpar[:,:,3:3+info[p]['nAntenna']])*np.exp(1.j * calpar[:,:,3+info[p]['nAntenna']:3+2*info[p]['nAntenna']] * math.pi / 180)).transpose((0,2,1))

		if os.path.isfile(additivefilenames[p]):
			adds.append(np.fromfile(additivefilenames[p], dtype='complex64').reshape((ttotal, nfreq, len(info[p]['subsetbl']))))
		else:
			adds.append(np.zeros((ttotal, nfreq, len(info[p]['subsetbl']))))

	#########start processing#######################
	t = []
	timing = []
	#datapulled = False
	for uvfile in uvfilenames:
		uvi = ap.miriad.UV(uvfile)
		if len(timing) > 0:
			print FILENAME + METHODNAME + "MSG:", uvfile + ' after', timing[-1]#uv.nchan

		if oppath == None:
			oppath = os.path.abspath(os.path.dirname(os.path.dirname(uvfile + '/'))) + '/'
		opuvname = oppath + os.path.basename(os.path.dirname(uvfile+'/')) + ano + 'O'
		print FILENAME + METHODNAME + "MSG: Creating %s"%opuvname
		if overwrite and os.path.isdir(opuvname):
			shutil.rmtree(opuvname)
		uvo = ap.miriad.UV(opuvname, status='new')
		uvo.init_from_uv(uvi)
		historystr = "Applied OMNICAL %s: "%time.asctime(time.localtime(time.time()))
		for cpfn, adfn in zip(calparfilenames, additivefilenames):
			historystr += os.path.abspath(cpfn) + ' ' + os.path.abspath(adfn) + ' '
		uvo['history'] += historystr + "\n"
		for preamble, data, flag in uvi.all(raw=True):
			uvo.copyvr(uvi)
			if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
				t += [preamble[1]]

				if len(t) > ttotal:
					print FILENAME + METHODNAME + " MSG: FATAL ERROR: calpar input array " + calparfilenames[p] + " has length", calpar.shape, "but the total length is exceeded when processing " + uvfile + " Aborted!"
					return
			polwanted = False
			for p, pol in zip(range(len(wantpols)), wantpols.keys()):
				if wantpols[pol] == uvi['pol']:
					a1, a2 = preamble[2]
					bl = bl1dmatrix[p][a1][a2]
					if bl > 0:
						additive = adds[p][len(t) - 1, :, bl - 1]
					elif bl < 0:
						additive = adds[p][len(t) - 1, :, - bl - 1].conjugate()
					else:
						additive = 0
						flag[:] = True
					#print data.shape, additive.shape, calpars[p][len(t) - 1, a1].shape
					uvo.write(preamble, (data-additive)/calpars[p][len(t) - 1, a1].conjugate()/calpars[p][len(t) - 1, a2], flag)
					polwanted = True
					break
			if not polwanted:
				uvo.write(preamble, data, flag)

		del(uvo)
		del(uvi)
		#if not datapulled:
			#print FILENAME + METHODNAME + " MSG:",  "FATAL ERROR: no data pulled from " + uvfile + ", check polarization information! Exiting."
			#exit(1)
	return


def stdmatrix(length, polydegree):#to find out the error in fitting y by a polynomial poly(x), one compute error vector by (I-A.(At.A)^-1 At).y, where Aij = i^j. This function returns (I-A.(At.A)^-1 At)
	A = np.array([[i**j for j in range(polydegree + 1)] for i in range(length)], dtype='int')
	At = A.transpose()
	return np.identity(length) - A.dot(la.pinv(At.dot(A), cond = 10**(-6)).dot(At))

#input two different redundant info, output True if they are the same and False if they are different
def compare_info(info1,info2, verbose=True, tolerance = 10**(-5)):
	try:
		floatkeys=['antloc','ubl','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
		intkeys = ['nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix']
		infomatrices=['A','B','At','Bt']
		specialkeys = ['ublindex']
		allkeys= floatkeys + intkeys + infomatrices + specialkeys#['antloc','ubl','nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB','A','B','At','Bt']
		diff=[]
		#10**5 for floating point errors
		for key in floatkeys:
			try:
				diff.append(round(la.norm(np.array(info1[key])-np.array(info2[key]))/tolerance)==0)
			except:
				diff.append(False)
		for key in intkeys:
			try:
				diff.append(la.norm(np.array(info1[key])-np.array(info2[key]))==0)
			except:
				diff.append(False)
		for key in infomatrices:
			try:
				diff.append(la.norm((info1[key]-info2[key]).todense())==0)
			except:
				diff.append(False)

		diff.append(True)
		try:
			for i,j in zip(info1['ublindex'],info2['ublindex']):
				diff[-1] = diff[-1] and (la.norm(np.array(i) - np.array(j))==0)
		except:
			diff[-1] = False
		bool = True
		for i in diff:
			bool = bool and i
		#print the first key found different (this will only trigger when the two info's have the same shape, so probably not very useful)
		if verbose and bool == False:
			for i in range(len(diff)):
				if diff[i] == False:
					print allkeys[i]
		return bool
	except ValueError:
		print "info doesn't have the same shape"
		return False

def omnical2omnigain(omnicalPath, utctimePath, info, outputPath = None):#outputPath should be a path without extensions like .omnigain which will be appended
	if outputPath == None:
		outputPath = omnicalPath.replace('.omnical', '')
	julDelta = 2415020
	#info = redundantCalibrator.info


	with open(utctimePath) as f:
		utctimes = f.readlines()
	calpars = np.fromfile(omnicalPath, dtype='float32')

	nT = len(utctimes)
	nF = len(calpars) / nT / (3 + 2 * info['nAntenna'] + 2 * info['nUBL'])
	#if nF != redundantCalibrator.nFrequency:
		#raise Exception('Error: time and frequency count implied in the infput files (%d %d) does not agree with those speficied in redundantCalibrator (%d %d). Exiting!'%(nT, nF, redundantCalibrator.nTime, redundantCalibrator.nFrequency))
	calpars = calpars.reshape((nT, nF, (3 + 2 * info['nAntenna'] + 2 * info['nUBL'])))

	jd = np.zeros((len(utctimes), 2), dtype='float32')#Julian dat is the only double in this whole thing so im storing it in two chunks as float
	sa = ephem.Observer()
	for utctime, t in zip(utctimes, range(len(utctimes))):
		sa.date = utctime
		jd[t, :] = struct.unpack('ff', struct.pack('d', sa.date + julDelta))

	opchisq = np.zeros((nT, 2 + 1 + 2 * nF), dtype = 'float32')
	opomnigain = np.zeros((nT, info['nAntenna'], 2 + 1 + 1 + 2 * nF), dtype = 'float32')
	opomnifit = np.zeros((nT, info['nUBL'], 2 + 3 + 1 + 2 * nF), dtype = 'float32')

	opchisq[:, :2] = jd
	opchisq[:, 2] = float(nF)
	opchisq[:, 3::2] = calpars[:, :, 0]#number of lincal iters
	opchisq[:, 4::2] = calpars[:, :, 2]#chisq which is sum of squares of errors in each visbility

	opomnigain[:, :, :2] = jd[:, None]
	opomnigain[:, :, 2] = np.array(info['subsetant']).astype('float32')
	opomnigain[:, :, 3] = float(nF)
	gains = (10**calpars[:, :, 3:(3 + info['nAntenna'])] * np.exp(1.j * math.pi * calpars[:, :, (3 + info['nAntenna']):(3 + 2 * info['nAntenna'])] / 180)).transpose((0,2,1))
	opomnigain[:, :, 4::2] = np.real(gains)
	opomnigain[:, :, 5::2] = np.imag(gains)

	opomnifit[:, :, :2] = jd[:, None]
	opomnifit[:, :, 2:5] = np.array(info['ubl']).astype('float32')
	opomnifit[:, :, 5] = float(nF)
	opomnifit[:, :, 6::2] = calpars[:, :, 3 + 2 * info['nAntenna']::2].transpose((0,2,1))
	opomnifit[:, :, 7::2] = calpars[:, :, 3 + 2 * info['nAntenna'] + 1::2].transpose((0,2,1))


	opchisq.tofile(outputPath + '.omnichisq')
	opomnigain.tofile(outputPath + '.omnigain')
	opomnifit.tofile(outputPath + '.omnifit')

class RedundantInfo(_O.RedundantInfo):#a class that contains redundant calibration information that should only be passed into C++
	def __init__(self, info):
		_O.RedundantInfo.__init__(self)
		if type(info) == type('a'):
			info = read_redundantinfo(info)
		elif type(info) != type({}):
			raise Exception("Error: info argument not recognized. It must be of either dictionary type (an info dictionary) *OR* string type (path to the info file).")
		
		for key in info.keys():
			try:
				if key in ['At','Bt']:
					tmp = []
					for i in range(info[key].shape[0]):
						for j in range(info[key].shape[1]):
							if info[key][i,j] != 0:
								tmp += [[i, j, info[key][i,j]]]
					self.__setattr__(key+'sparse', np.array(tmp, dtype = 'int32'))
				elif key in ['A','B']:
					self.__setattr__(key, info[key].todense().astype('int32'))
				elif key in ['ublindex']:
					tmp = []
					for i in range(len(info[key])):
						for j in range(len(info[key][i])):
							tmp += [[i, j, info[key][i][j][0], info[key][i][j][1], info[key][i][j][2]]]
					self.__setattr__(key, np.array(tmp, dtype='int32'))
				elif key in int_infokeys:
					self.__setattr__(key, int(info[key]))
				elif key in intarray_infokeys and key != 'ublindex':
					self.__setattr__(key, np.array(info[key]).astype('int32'))
				elif key in float_infokeys:
					self.__setattr__(key, np.array(info[key]).astype('float32'))
			except:
				raise Exception("Error parsing %s item."%key)

	
	def get_info(self):
		info = {}
		for key in infokeys:
			try:
				if key in ['A','B']:
					#print key
					info[key] = sps.csr_matrix(self.__getattribute__(key))
				elif key in ['At','Bt']:
					tmp = self.__getattribute__(key+'sparse')
					matrix = np.zeros((info['nAntenna'] + info['nUBL'], len(info['crossindex'])))
					for i in tmp:
						matrix[i[0],i[1]] = i[2]
					info[key] = sps.csr_matrix(matrix)
				elif key in ['ublindex']:
					ublindex = []
					for i in self.__getattribute__(key):
						while len(ublindex) < i[0] + 1:
							ublindex.append(np.zeros((1,3)))
						while len(ublindex[i[0]]) < i[1] + 1:
							ublindex[i[0]] = np.array(ublindex[i[0]].tolist() + [[0,0,0]])
						ublindex[i[0]][i[1]] = np.array(i[2:])
					info[key] = ublindex

				else:
					#print key
					info[key] = self.__getattribute__(key)
			except:
				raise Exception("Error retrieving %s item."%key)
		return info
		
		
class RedundantCalibrator:
#This class is the main tool for performing redundant calibration on data sets. For a given redundant configuration, say 32 antennas with 3 bad antennas, the user should create one instance of Redundant calibrator and reuse it for all data collected from that array. In general, upon creating an instance, the user need to create the info field of the instance by either computing it or reading it from a text file. readyForCpp(verbose = True) should be a very helpful function to provide information on what information is missing for running the calibration.
	def __init__(self, nTotalAnt, info = None):
		methodName = '.__init__.'
		self.className = '.RedundantCalibrator.'
		self.nTotalAnt = nTotalAnt
		self.nTotalBaselineAuto = (self.nTotalAnt + 1) * self.nTotalAnt / 2
		self.nTotalBaselineCross = (self.nTotalAnt - 1) * self.nTotalAnt / 2
		self.antennaLocation = np.zeros((self.nTotalAnt, 3))
		self.antennaLocationTolerance = 10**(-6)
		self.badAntenna = []
		self.badUBL = []
		self.totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(self.nTotalAnt)])#PAPER miriad convention by default
		
		self.Info = None
		self.infoPath = './tmp_redundantinfo'
		self.infoFileExist = False
		self.dataFileExist = False
		self.keepData = False
		self.tmpDataPath = './tmp_calibration_omni_data'
		self.dataPath = self.tmpDataPath #complex128 type binary visibility file
		self.keepCalpar = False
		self.calparPath = None
		self.nFrequency = -1
		self.nTime = -1
		self.removeDegeneracy = True
		self.removeAdditive = False
		self.removeAdditivePeriod = -1
		self.convergePercent = 0.01 #convergence criterion in relative change of chi^2. By default it stops when reaches 0.01, namely 1% decrease in chi^2.
		self.maxIteration = 50 #max number of iterations in lincal
		self.stepSize = 0.3 #step size for lincal. (0, 1]. < 0.4 recommended.
		self.computeUBLFit = True
		
		if info != None:
			if type(info) == type({}):
				
				self.Info = RedundantInfo(info)
			elif type(info) == type('a'):
				self.read_redundantinfo(info)
			else:
				raise Exception(self.className + methodName + "Error: info argument not recognized. It must be of either dictionary type (an info dictionary) *OR* string type (path to the info file).")

	def read_redundantinfo(self, infopath):#redundantinfo is necessary for running redundant calibration. The text file should contain 29 lines each describes one item in the info.
		self.infoPath = infopath
		
		self.Info = RedundantInfo(read_redundantinfo(infopath))
		self.infoFileExist = True

	def write_redundantinfo(self, infoPath = None, overwrite = False):
		methodName = '.write_redundantinfo.'
		if infoPath == None:
			infoPath = self.infoPath
		if (self.Info != None) and (infoPath != None):
			write_redundantinfo(self.Info.get_info(), infoPath, overwrite = overwrite)
			self.infoPath = infoPath
			self.infoFileExist = True
		else:
			raise Exception(self.className + methodName + "Error: either 1) info does not yet exist for the current instance, or 2) an info file already exists on disk, or 3) no file path is ever specified.")

	def read_arrayinfo(self, arrayinfopath, verbose = False):#array info is the minimum set of information to uniquely describe a redundant array, and is needed to compute redundant info. It includes, in each line, bad antenna indices, bad unique baseline indices, tolerance of error when checking redundancy, antenna locations, and visibility's antenna pairing conventions. Unlike redundant info which is a self-contained dictionary, items in array info each have their own fields in the instance.
		methodName = ".read_arrayinfo."
		if not os.path.isfile(arrayinfopath):
			raise Exception(self.className + methodName + "Error: Array info file " + arrayinfopath + " doesn't exist!")
		with open(arrayinfopath) as f:
			rawinfo = [[float(x) for x in line.split()] for line in f]
		if verbose:
			print self.className + methodName + " MSG:",  "Reading", arrayinfopath, "...",

		self.badAntenna = np.array(rawinfo[0]).astype(int)
		if self.badAntenna[0] < 0:
			self.badAntenna = np.zeros(0)

		self.badUBL = np.array(rawinfo[1]).astype(int)
		if self.badUBL[0] < 0:
			self.badUBL = np.zeros(0)

		self.antennaLocationTolerance = rawinfo[2][0]

		for a in range(len(self.antennaLocation)):
			if len(rawinfo[a + 3]) != 3:
				raise Exception(self.className + methodName + "Error: Format error in " + arrayinfopath + ": The antenna locations should start on the 4th line, with 3 numbers in each line!")
			else:
				self.antennaLocation[a] = np.array(rawinfo[a + 3])

		for bl in range(len(self.totalVisibilityId)):
			if len(rawinfo[bl + 3 + len(self.antennaLocation)]) != 2:
				raise Exception(self.className + methodName + "Error: Format error in " + arrayinfopath + ": The baseline to antenna mapping should start after antenna locations, with 2 numbers (conj index, index) in each line!")
			else:
				self.totalVisibilityId[bl] = np.array(rawinfo[bl + 3 + len(self.antennaLocation)]).astype(int)
		if verbose:
			print "Bad antenna indices:", self.badAntenna,
			print "Bad UBL indices:", self.badUBL


	def readyForCpp(self, verbose = True):#check if all necessary parameters are specified to call Cpp
		methodName = '.readyForCpp.'
		#if not (self.dataFileExist and os.path.isfile(self.dataPath)):
			#if verbose:
				#print self.className + methodName + "Error: data file check failed. No data file specified or specified filename does not exist."
			#return False
		#if os.path.getsize(self.dataPath) / 8 != self.nTime * self.nFrequency * self.nTotalBaselineAuto:
			#if verbose:
				#print self.className + methodName + "Error: data size check failed. File on disk seems to contain " + str(os.path.getsize(self.dataPath) / 8) + " complex64 numbers, where as we expect " + str(self.nTime * self.nFrequency * self.nTotalBaselineAuto) + '.'
			#return False

		#if not self.infoFileExist :
			#if verbose:
				#print self.className + methodName + "Error: info file existence check failed. Call read_redundantinfo(self, infoPath) function to read in an existing redundant info text file or compute_redundantinfo(arrayinfoPath=None) to compute redundantinfo write_redundantinfo() to write redundantinfo."
			#return False
		if self.Info == None :
			if verbose:
				print self.className + methodName + "Error: self.Info existence check failed."
			return False

		#if self.removeAdditive and self.removeAdditivePeriod <= 0:
			#if verbose:
				#print self.className + methodName + "Error: removeAdditive option is True but the removeAdditivePeriod parameter is negative (invalid)."
			#return False

		if abs(self.convergePercent - 0.5) >= 0.5 or self.maxIteration <= 0 or abs(self.stepSize - 0.5) >= 0.5:
			if verbose:
				print self.className + methodName + "Error: lincal parameter check failed. convergePercent and stepSize should be between 0 and 1, and maxIteration has to be positive integer."
			return False
		if verbose:
			print self.className + methodName + "Check passed."
		return True

	def cal(self, data, additivein, verbose = False):#data can be either 3d numpy array or a binary file path
		methodName = '.cal.'
		if type(data) == type(' '):
			raise Exception("passing in file on disk no longer supported!")
		elif type(data) == type(np.zeros(1)) and len(data.shape) == 3 and len(data[0,0]) == self.Info.nBaseline and self.dataPath == self.tmpDataPath:
			(self.nTime, self.nFrequency, _) = data.shape
			np.array(data, dtype = 'complex64').tofile(self.tmpDataPath)
			self.dataPath = self.tmpDataPath
			self.dataFileExist = True
		else:
			raise Exception(self.className + methodName + "Error: data type must be a file path name to a binary file *OR* a 3D numpy array of dimensions (nTime, nFrequency, nTotalBaselineAuto). You have either 1) passed in the wrong type, or 2) passed in a correct data array but have mismatching self.dataPath and self.tmpDataPath (these paths will be overwritten if you pass in an array as data!).")

		#if (not self.infoFileExist):
			#self.write_redundantinfo()

		if self.readyForCpp(verbose = False):
			#print data.shape, self.rawCalpar.shape, self.Info.nUBL, additivein.shape, int(self.removeDegeneracy), int(self.calMode), float(self.convergePercent), int(self.maxIteration), float(self.stepSize)
			_ =_O.redcal(data, self.rawCalpar, self.Info, additivein, removedegen = int(self.removeDegeneracy), uselogcal = int(self.calMode),maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize))
			#, removedegen = int(self.removeDegeneracy), uselogcal = int(self.calMode), tol = float(self.convergePercent), maxiter=int(self.maxIteration), stepsize=float(self.stepSize)
			#_O.omnical(self.dataPath, self.infoPath, int(self.nTime), int(self.nFrequency), int(self.nTotalAnt), int(self.removeDegeneracy), int(self.removeAdditive), str(self.removeAdditivePeriod), int(self.calMode), float(self.convergePercent), int(self.maxIteration), float(self.stepSize))

			#command = "./omnical " + self.dataPath + " " + self.infoPath + " " + str(self.nTime) + " " + str(self.nFrequency) + " "  + str(self.nTotalAnt) + " " + str(int(self.removeDegeneracy)) + " " + str(int(self.removeAdditive)) + " " + str(self.removeAdditivePeriod) + " " + self.calMode + " " + str(self.convergePercent) + " " + str(self.maxIteration) + " " + str(self.stepSize)
			#if verbose:
				#print self.className + methodName + "System call: " + command
			#os.system(command)

			if self.removeAdditive and self.removeAdditivePeriod > 0:
				self.calparPath = self.dataPath + '_add' + str(self.removeAdditivePeriod) + '.omnical'
			else:
				self.calparPath = self.dataPath + '.omnical'
			#self.rawCalpar = np.fromfile(self.calparPath, dtype = 'float32').reshape((self.nTime, self.nFrequency, 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)))
			if self.calMode == '0' or self.calMode == '1':
				self.chisq = self.rawCalpar[:, :, 2]
			elif self.calMode == '2':
				self.chisq = self.rawCalpar[:, :, 1]
			self.calpar = np.zeros((self.nTime, self.nFrequency, self.nTotalAnt), dtype='complex64')
			self.calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
			self.bestfit = self.rawCalpar[:, :, (3 + 2 * self.Info.nAntenna):: 2] + 1.j * self.rawCalpar[:, :, (4 + 2 * self.Info.nAntenna):: 2]
			#if not self.keepCalpar:
				#os.remove(self.calparPath)
			#if not self.keepData and self.dataPath == self.tmpDataPath:
				#os.remove(self.dataPath)
		else:
			raise Exception(self.className + methodName + "Error: function is called prematurely. The current instance failed the readyForCpp() check. Try instance.readyForCpp() for more info.")


	def loglincal(self, data, additivein, verbose = False):
		self.calMode = '1'
		self.cal(data, additivein, verbose)

	def lincal(self, data, verbose = False):
		if self.rawCalpar.shape != (len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)):
			raise Exception("ERROR: lincal called without a properly shaped self.rawCalpar! Excpeted shape is (%i, %i, %i)!"%(len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)))
		_O.redcal(data, self.rawCalpar, self.Info, additivein, removedegen = int(self.removeDegeneracy), uselogcal = 0, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit))
		self.chisq = self.rawCalpar[:, :, 2]
		self.calpar = np.zeros((self.nTime, self.nFrequency, self.nTotalAnt), dtype='complex64')
		self.calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
		self.bestfit = self.rawCalpar[:, :, (3 + 2 * self.Info.nAntenna):: 2] + 1.j * self.rawCalpar[:, :, (4 + 2 * self.Info.nAntenna):: 2]
			
	def logcal(self, data, verbose = False):
		self.rawCalpar = np.zeros((len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)), dtype = 'float32')
		_O.redcal(data, self.rawCalpar, self.Info, additivein, removedegen = int(self.removeDegeneracy), uselogcal = 1, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit))
		self.chisq = self.rawCalpar[:, :, 2]
		self.calpar = np.zeros((self.nTime, self.nFrequency, self.nTotalAnt), dtype='complex64')
		self.calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
		self.bestfit = self.rawCalpar[:, :, (3 + 2 * self.Info.nAntenna):: 2] + 1.j * self.rawCalpar[:, :, (4 + 2 * self.Info.nAntenna):: 2]

	def set_badUBL(self, badUBL):
		if np.array(badUBL).shape[-1] != 3 or len(np.array(badUBL).shape) != 2:
			raise Exception("ERROR: badUBL need to be a list of coordinates each with 3 numbers.")
		badindex = []
		UBL = self.compute_UBL(self.antennaLocationTolerance)
		for badubl in badUBL:
			for i, ubl in zip(range(len(UBL)), UBL):
				if la.norm(badubl - ubl) < self.antennaLocationTolerance:
					badindex += [i]
					break
		if len(badindex) != len(badUBL):
			raise Exception("ERROR: some badUBL not found in self.computeUBL!")
		else:
			self.badUBL = np.sort(badindex).tolist()

	def compute_redundantinfo(self, arrayinfoPath = None):
		if arrayinfoPath != None and os.path.isfile(arrayinfoPath):
			self.read_arrayinfo(arrayinfoPath)
		if np.linalg.norm(self.antennaLocation) == 0:
			raise Exception("Error: compute_redundantinfo() called before self.antennaLocation is specified. Use configFilePath option when calling compute_redundantinfo() to specify array info file, or manually set self.antennaLocation for the RedundantCalibrator instance.")

		#nAntenna and subsetant : get rid of the bad antennas
		nant=len(self.antennaLocation)
		subsetant=[i for i in range(nant) if i not in self.badAntenna]
		nAntenna=len(subsetant)
		antloc=[self.antennaLocation[ant] for ant in subsetant]
		##########################################################################################
		#find out ubl
		#use the function compute_UBL to find the ubl
		tolerance=self.antennaLocationTolerance;
		ublall=self.compute_UBL(tolerance)
		#delete the bad ubl's
		ubl=np.delete(ublall,np.array(self.badUBL),0)
		nUBL=len(ubl);
		#################################################################################################
		#calculate the norm of the difference of two vectors (just la.norm actually)
		def dis(a1,a2):
			return np.linalg.norm(np.array(a1)-np.array(a2))
		#find nBaseline (include auto baselines) and subsetbl
		badbl=[ublall[i] for i in self.badUBL]
		nbl=0;
		goodpairs=[];
		for i in range(len(antloc)):
			for j in range(i+1):
				bool=False
				for bl in badbl:
					bool = bool or dis(antloc[i]-antloc[j],bl)<tolerance or dis(antloc[i]-antloc[j],-bl)<tolerance
				if bool == False:
					nbl+=1
					goodpairs.append([i,j])
		nBaseline=len(goodpairs)
		#from a pair of good antenna index to baseline index
		subsetbl=np.array([self.get_baseline([subsetant[bl[0]],subsetant[bl[1]]]) for bl in goodpairs])
		##################################################################################
		#bltoubl: cross bl number to ubl index
		def findublindex(pair,ubl=ubl):
			i=pair[0]
			j=pair[1]
			for k in range(len(ubl)):
				if dis(antloc[i]-antloc[j],ubl[k])<tolerance or dis(antloc[i]-antloc[j],-ubl[k])<tolerance:
					return k
			return "no match"
		bltoubl=[];
		for i in goodpairs:
			if i[0]!=i[1]:
				bltoubl.append(findublindex(i))
		#################################################################################
		#reversed:   cross only bl if reversed -1, otherwise 1
		crosspair=[]
		for p in goodpairs:
			if p[0]!=p[1]:
				crosspair.append(p)
		reverse=[]
		for k in range(len(crosspair)):
			i=crosspair[k][0]
			j=crosspair[k][1]
			if dis(antloc[i]-antloc[j],ubl[bltoubl[k]])<tolerance:
				reverse.append(1)
			elif dis(antloc[i]-antloc[j],-ubl[bltoubl[k]])<tolerance:
				reverse.append(-1)
			else :
				print "something's wrong with bltoubl"
		######################################################################################
		#reversedauto: the index of good baselines (auto included) in all baselines
		#autoindex: index of auto bls among good bls
		#crossindex: index of cross bls among good bls
		#ncross
		reversedauto = range(len(goodpairs))
		#find the autoindex and crossindex in goodpairs
		autoindex=[]
		crossindex=[]
		for i in range(len(goodpairs)):
			if goodpairs[i][0]==goodpairs[i][1]:
				autoindex.append(i)
			else:
				crossindex.append(i)
		for i in autoindex:
			reversedauto[i]=1
		for i in range(len(crossindex)):
			reversedauto[crossindex[i]]=reverse[i]
		reversedauto=np.array(reversedauto)
		autoindex=np.array(autoindex)
		crossindex=np.array(crossindex)
		ncross=len(crossindex)
		###################################################
		#bl2d:  from 1d bl index to a pair of antenna numbers
		bl2d=[]
		for pair in goodpairs:
			bl2d.append(pair[::-1])
		bl2d=np.array(bl2d)
		###################################################
		#ublcount:  for each ubl, the number of good cross bls corresponding to it
		countdict={}
		for bl in bltoubl:
			countdict[bl]=0

		for bl in bltoubl:
			countdict[bl]+=1

		ublcount=[]
		for i in range(nUBL):
			ublcount.append(countdict[i])
		ublcount=np.array(ublcount)
		####################################################################################
		#ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
		countdict={}
		for bl in bltoubl:
			countdict[bl]=[]

		for i in range(len(crosspair)):
			ant1=crosspair[i][1]
			ant2=crosspair[i][0]
			countdict[bltoubl[i]].append([ant1,ant2,i])

		ublindex=[]
		for i in range(nUBL):
			ublindex.append(countdict[i])
		#turn each list in ublindex into np array
		for i in range(len(ublindex)):
			ublindex[i]=np.array(ublindex[i])
		ublindex=np.array(ublindex)
		###############################################################################
		#bl1dmatrix: a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
				#I suppose 99999 for bad and auto baselines?
		bl1dmatrix=99999*np.ones([nAntenna,nAntenna],dtype='int16')
		for i in range(len(crosspair)):
			bl1dmatrix[crosspair[i][1]][crosspair[i][0]]=i
			bl1dmatrix[crosspair[i][0]][crosspair[i][1]]=i
		####################################################################################3
		#degenM:
		a=[]
		for i in range(len(antloc)):
			a.append(np.append(antloc[i],1))
		a=np.array(a)

		d=[]
		for i in range(len(ubl)):
			d.append(np.append(ubl[i],0))
		d=np.array(d)

		m1=-a.dot(la.pinv(np.transpose(a).dot(a), cond = 10**(-6))).dot(np.transpose(a))
		m2=d.dot(la.pinv(np.transpose(a).dot(a), cond = 10**(-6))).dot(np.transpose(a))
		degenM = np.append(m1,m2,axis=0)
		#####################################################################################
		#A: A matrix for logcal amplitude
		A=np.zeros([len(crosspair),nAntenna+len(ubl)])
		for i in range(len(crosspair)):
			A[i][crosspair[i][0]]=1
			A[i][crosspair[i][1]]=1
			A[i][nAntenna+bltoubl[i]]=1
		A=sps.csr_matrix(A)
		#################################################################################
		#B: B matrix for logcal phase
		B=np.zeros([len(crosspair),nAntenna+len(ubl)])
		for i in range(len(crosspair)):
			B[i][crosspair[i][0]]=reverse[i]*1
			B[i][crosspair[i][1]]=reverse[i]*-1
			B[i][nAntenna+bltoubl[i]]=1
		B=sps.csr_matrix(B)
		###########################################################################
		#create info dictionary
		info={}
		info['nAntenna']=nAntenna
		info['nUBL']=nUBL
		info['nBaseline']=nBaseline
		info['subsetant']=subsetant
		info['antloc']=antloc
		info['subsetbl']=subsetbl
		info['ubl']=ubl
		info['bltoubl']=bltoubl
		info['reversed']=reverse
		info['reversedauto']=reversedauto
		info['autoindex']=autoindex
		info['crossindex']=crossindex
		#info['ncross']=ncross
		info['bl2d']=bl2d
		info['ublcount']=ublcount
		info['ublindex']=ublindex
		info['bl1dmatrix']=bl1dmatrix
		info['degenM']=degenM
		info['A']=A
		info['B']=B
		with warnings.catch_warnings():
				warnings.filterwarnings("ignore",category=DeprecationWarning)
				info['At'] = info['A'].transpose()
				info['Bt'] = info['B'].transpose()
				info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense(), cond = 10**(-6))#(AtA)^-1
				info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense(), cond = 10**(-6))#(BtB)^-1
				info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
				info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
				info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
				info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
				info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
				info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB
		self.Info = RedundantInfo(info)



	#inverse function of totalVisibilityId, calculate the baseline index from the antenna pair
	def get_baseline(self,pair):
		if not (type(pair) == list or type(pair) == np.ndarray or type(pair) == tuple):
			raise Exception("input needs to be a list of two numbers")
			return
		elif len(np.array(pair)) != 2:
			raise Exception("input needs to be a list of two numbers")
			return
		elif type(pair[0]) == str or type(pair[0]) == np.string_:
			raise Exception("input needs to be number not string")
			return
		sortp = np.array(sorted(pair))
		for i in range(len(self.totalVisibilityId)):
			if self.totalVisibilityId[i][0] == sortp[0] and self.totalVisibilityId[i][1] == sortp[1]:
				return i
		raise Exception("antenna index out of range")

	#with antenna locations and tolerance, calculate the unique baselines. (In the order of omniscope baseline index convention)
	#compute_UBL returns the average of all baselines in that ubl group
	def compute_UBL(self,tolerance = 0.1):
		#check if the tolerance is not a string
		if type(tolerance) == str:
			raise Exception("tolerance needs to be number not string")
			return
		#remove the bad antennas
		nant=len(self.antennaLocation)
		subsetant=[i for i in range(nant) if i not in self.badAntenna]
		nAntenna=len(subsetant)
		antloc = np.array([self.antennaLocation[ant] for ant in subsetant])
		ubllist = np.array([np.array([np.array([0,0,0]),1])]);
		for i in range(len(antloc)):
			for j in range(i+1,len(antloc)):
				bool = True
				for k in range(len(ubllist)):
					if  la.norm(antloc[i]-antloc[j]-ubllist[k][0])<tolerance:
						n=ubllist[k][1]
						ubllist[k][0]=1/(n+1.0)*(n*ubllist[k][0]+antloc[i]-antloc[j])
						ubllist[k][1]+=1
						bool = False
					elif  la.norm(antloc[i]-antloc[j]+ubllist[k][0])<tolerance:
						n=ubllist[k][1]
						ubllist[k][0]=1/(n+1.0)*(n*ubllist[k][0]-(antloc[i]-antloc[j]))
						ubllist[k][1]+=1
						bool = False
				if bool :
					ubllist = np.append(ubllist,np.array([np.array([antloc[j]-antloc[i],1])]),axis=0)
		ubllist = np.delete(ubllist,0,0)
		ublall=[]
		for ubl in ubllist:
			ublall.append(ubl[0])
		ublall=np.array(ublall)
		return ublall

	#need to do compute_redundantinfo first for this function to work (needs 'bl1dmatrix')
	#input the antenna pair(as a list of two numbers), return the corresponding ubl index
	def get_ublindex(self,antpair):
		#check if the input is a list, tuple, np.array of two numbers
		if not (type(antpair) == list or type(antpair) == np.ndarray or type(antpair) == tuple):
			raise Exception("input needs to be a list of two numbers")
			return
		elif len(np.array(antpair)) != 2:
			raise Exception("input needs to be a list of two numbers")
			return
		elif type(antpair[0]) == str or type(antpair[0]) == np.string_:
			raise Exception("input needs to be number not string")
			return

		#check if self.info['bl1dmatrix'] exists
		try:
			_ = self.Info.bl1dmatrix
		except:
			raise Exception("needs Info.bl1dmatrix")

		crossblindex=self.Info.bl1dmatrix[antpair[0]][antpair[1]]
		if antpair[0]==antpair[1]:
			return "auto correlation"
		elif crossblindex == 99999:
			return "bad ubl"
		return self.Info.bltoubl[crossblindex]


	#need to do compute_redundantinfo first
	#input the antenna pair, return -1 if it is a reversed baseline and 1 if it is not reversed
	def get_reversed(self,antpair):
		#check if the input is a list, tuple, np.array of two numbers
		if not (type(antpair) == list or type(antpair) == np.ndarray or type(antpair) == tuple):
			raise Exception("input needs to be a list of two numbers")
			return
		elif len(np.array(antpair)) != 2:
			raise Exception("input needs to be a list of two numbers")
			return
		elif type(antpair[0]) == str or type(antpair[0]) == np.string_:
			raise Exception("input needs to be number not string")
			return

		#check if self.info['bl1dmatrix'] exists
		try:
			_ = self.Info.bl1dmatrix
		except:
			raise Exception("needs Info.bl1dmatrix")

		crossblindex=self.Info.bl1dmatrix[antpair[0]][antpair[1]]
		if antpair[0] == antpair[1]:
			return 1
		if crossblindex == 99999:
			return 'badbaseline'
		return self.Info.reversed[crossblindex]







