import datetime
import socket, multiprocessing, math, random, traceback, ephem, string, commands, datetime
import time
from time import ctime
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

def read_redundantinfo(infopath):
	with open(infopath) as f:
		rawinfo = np.array([np.array([float(x) for x in line.split()]) for line in f])
	METHODNAME = "read_redundantinfo"
	info = {}
	infocount = 0;
	info['nAntenna'] = int(rawinfo[infocount][0]) #number of good antennas among all (64) antennas, same as the length of subsetant
	infocount += 1

	info['nUBL'] = int(rawinfo[infocount][0]) #number of unique baselines
	infocount += 1
		
	nbl = int(rawinfo[infocount][0])
	info['nbl'] = nbl
	infocount += 1
	ncross = nbl - info['nAntenna']
	info['ncross'] = ncross

	info['subsetant'] = rawinfo[infocount].astype(int) #the index of good antennas in all (64) antennas
	infocount += 1
	info['subsetbl'] = rawinfo[infocount].astype(int) #the index of good baselines (auto included) in all baselines
	infocount += 1
	info['ubl'] = rawinfo[infocount].reshape((info['nUBL'],3)).astype(int) #unique baseline vectors
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

	return info

antlocvecX5 = np.array([[0.999973, -0.0105871, 0.00142779], [0.010705, 
  0.999799, -0.00475158], [-0.00263492, 0.0167222, -0.0178542]], dtype='float32')

def stdmatrix(length, polydegree):#to find out the error in fitting y by a polynomial poly(x), one compute error vector by (I-A.(At.A)^-1 At).y, where Aij = i^j. This function returns (I-A.(At.A)^-1 At)
	A = np.array([[i**j for j in range(polydegree + 1)] for i in range(length)], dtype='int')
	At = A.transpose()
	return np.identity(length) - A.dot(la.pinv(At.dot(A)).dot(At))
