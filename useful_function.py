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
    
import calibration_omni as omni

calibrator = omni.RedundantCalibrator(32)

#inverse function of totalVisibilityId, calculate the baseline index from the antenna pair
def toBaseline(self,pair):
	sortp=np.array(sorted(pair))
	for i in range(len(tvid)):
		if self.totalVisibilityId[i][0] == sortp[0] and self.totalVisibilityId[i][1] == sortp[1]:
			return i
	return 'no match'

#with antenna locations and tolerance, calculate the unique baselines. (In the order of omniscope baseline index convention)
def UBL(self,tolerance):
	antloc=self.antennaLocation
	ubllist=np.array([np.array([0,0,0])]);
	for i in range(len(antloc)):
		for j in range(i+1,len(antloc)):
			bool = False;
			for bl in ubllist:
				bool = bool or (la.norm(antloc[i]-antloc[j]-bl)<tolerance or la.norm(antloc[i]-antloc[j]+bl)<tolerance)
			if bool == False:			
				ubllist = np.concatenate((ubllist,[antloc[j]-antloc[i]]))
	ublall = np.delete(ubllist,0,0)
	return ublall

#need to do compute_info first for this function to work
#input the antenna pair, return the corresponding ubl index
def pair2ublindex(self,antpair):
	crossblindex=self.info['bl1dmatrix'][antpair[0]][antpair[1]]
	if antpair[0]==antpair[1]:
		return "auto correlation"
	elif crossblindex == 99999:
		return "bad ubl"
	return self.info['bltoubl'][crossblindex]
	
def pair2
	
	
	
	
