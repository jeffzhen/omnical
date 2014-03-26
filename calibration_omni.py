import datetime
import socket, multiprocessing, math, random, traceback, ephem, string, commands, datetime
import time
from time import ctime
import aipy as ap
import struct
import numpy as np
import os, sys
import datetime
from optparse import OptionParser
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy as sp
    import scipy.sparse as sps
    import scipy.linalg as la

FILENAME = "calibration_omni.py"

def read_redundantinfo(infopath):
	with open(infopath) as f:
		rawinfo = np.array([np.array([float(x) for x in line.split()]) for line in f])
	METHODNAME = "read_redundantinfo"
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
	
	
	info['bl1dmatrix'] = rawinfo[infocount].reshape((info['nAntenna'], info['nAntenna'])).astype(int) #a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index counting auto corr
	infocount += 1
	#matrices
	info['degenM'] = rawinfo[infocount].reshape((info['nAntenna'] + info['nUBL'], info['nAntenna']))
	infocount += 1	
	info['A'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #A matrix for logcal amplitude
	infocount += 1
	info['B'] = sps.csr_matrix(rawinfo[infocount].reshape((ncross, info['nAntenna'] + info['nUBL'])).astype(int)) #B matrix for logcal amplitude
	infocount += 1
	##The sparse matrices are treated a little differently because they are not rectangular
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore",category=DeprecationWarning)
		info['At'] = info['A'].transpose()
		info['Bt'] = info['B'].transpose()
		info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense())#(AtA)^-1
		info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense())#(BtB)^-1
		info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
		info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
		info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
		info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
		info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
		info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB
	print "done. nAntenna, nUBL, nBaseline = ", len(info['subsetant']), info['nUBL'], info['nBaseline']
	return info

antlocvecX5 = np.array([[0.999973, -0.0105871, 0.00142779], [0.010705, 
  0.999799, -0.00475158], [-0.00263492, 0.0167222, -0.0178542]], dtype='float32')

def importuvs(uvfilenames, info, wantpols):
	METHODNAME = "*importuvs*"
	latP = -0.53619181096511903#todo: figure this out from uv files
	lonP = 0.37399448506783717
	############################################################
	sa = ephem.Observer()
	sa.lon = lonP #//todo: read from uv property
	sa.lat = latP
	sun = ephem.Sun()
	julDelta = 2415020 # =julian date - pyephem's Observer date
	####get some info from the first uvfile
	uv=ap.miriad.UV(uvfilenames[0])
	nfreq = uv.nchan;
	nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
	startfreq = uv['sfreq']
	dfreq = uv['sdf']
	del(uv)
	####prepare processing
	deftime = 2000
	data = np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
	#sunpos = np.zeros((deftime, 2))
	t = []
	timing = []
	lst = []
	
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
				lst += [(float(sa.sidereal_time()) * 24./2./math.pi)]
				if len(t) > len(data):
					print FILENAME + METHODNAME + " MSG:",  "expanding number of time slices from", len(data), "to", len(data) + deftime
					data = np.concatenate((data, np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64'))) 
					#sunpos = np.concatenate((sunpos, np.zeros((deftime, 2))))
					#sunpos[len(t) - 1] = np.asarray([[sun.alt, sun.az]])
			for p, pol in zip(range(len(wantpols)), wantpols.keys()):
				if wantpols[pol] == uv['pol']:#//todo: use select() 
					a1, a2 = preamble[2]
					bl = info[p]['bl1dmatrix'][a1, a2]
					if bl < info[p]['nBaseline']:
						datapulled = True
						#print info[p]['subsetbl'][info[p]['crossindex'][bl]],
						data[len(t) - 1, p, info[p]['subsetbl'][info[p]['crossindex'][bl]]] = rawd.data.astype('complex64')
		del(uv)
		if not datapulled:
			print FILENAME + METHODNAME + " MSG:",  "FATAL ERROR: no data pulled from " + uvfile + ", check polarization information! Exiting."
			exit(1)
	return data[:len(t)], t, timing, lst

def apply_omnical_uvs(uvfilenames, calparfilenames, info, wantpols, oppath, ano):
	METHODNAME = "*apply_omnical_uvs*"

	####get some info from the first uvfile
	uv=ap.miriad.UV(uvfilenames[0])
	nfreq = uv.nchan;
	nant = uv['nants'] / 2 # 'nants' counting ant-pols, so divide 2
	startfreq = uv['sfreq']
	dfreq = uv['sdf']
	del(uv)
	
	####load calpar and check dimensions, massage calpar from txfx(3+2a+2u) to t*goodabl*f
	blcalpar = []#calpar for each baseline, auto included
	for p in range(len(wantpols)):
		calpar = np.fromfile(calparfilenames[p], dtype='float32')	
		if len(calpar)%(nfreq *( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']))) != 0:
			print FILENAME + METHODNAME + " MSG:",  "FATAL ERROR: calpar input array " + calparfilenames[p] + " has length", calpar.shape, "which is not divisible by ", nfreq, 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']), "Aborted!"
			return
		ttotal = len(calpar)/(nfreq *( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL'])))
		calpar = calpar.reshape((ttotal, nfreq, ( 3 + 2 * (info[p]['nAntenna'] + info[p]['nUBL']))))
		calpar = (10**calpar[:,:,3:3+info[p]['nAntenna']])*np.exp(1.j*calpar[:,:,3+info[p]['nAntenna']:3+2*info[p]['nAntenna']])
		blcalpar.append(1 + np.zeros((ttotal, info[p]['nBaseline'], nfreq),dtype='complex64'))
		for bl in range(info[p]['nBaseline']):
			blcalpar[p][:, bl, :] *= (calpar[:, :, info[p]['bl2d'][bl,0]].conj() * calpar[:, :, info[p]['bl2d'][bl, 1]])



	#########start processing#######################
	t = []
	timing = []
	#datapulled = False
	for uvfile in uvfilenames:
		uvi = ap.miriad.UV(uvfile)
		if len(timing) > 0:	
			print FILENAME + METHODNAME + "MSG:", uvfile + ' after', timing[-1]#uv.nchan
		uvo = ap.miriad.UV(oppath + os.path.basename(os.path.dirname(uvfile+'/')) + ano + 'omnical', status='new')
		uvo.init_from_uv(uvi)
		historystr = "Applied "
		for cpfn in calparfilenames:
			historystr += cpfn
#		uvo.pipe(uvi, mfunc=applycp, append2hist=historystr + "\n")
		for preamble, data, flag in uvi.all(raw=True):
			uvo.copyvr(uvi)
			if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
				t += [preamble[1]]

				if len(t) > ttotal:
					print FILENAME + METHODNAME + " MSG: FATAL ERROR: calpar input array " + calparfilenames[p] + " has length", calpar.shape, "but the total length is exceeded when processing " + uvfile + " Aborted!"
					return
			for p, pol in zip(range(len(wantpols)), wantpols.keys()):
				if wantpols[pol] == uvi['pol']:
					a1, a2 = preamble[2]
					bl = info[p]['bl1dmatrix'][a1, a2]
					if bl < info[p]['ncross']:
						#datapulled = True
						#print info[p]['subsetbl'][info[p]['crossindex'][bl]],
						uvo.write(preamble, data/blcalpar[p][len(t) - 1, info[p]['crossindex'][bl]], flag)
					#//todo: correct autocorr as well
					else:
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
	return np.identity(length) - A.dot(la.pinv(At.dot(A)).dot(At))
