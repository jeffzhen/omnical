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

correctinfo = omni.read_redundantinfo('redundantinfo_PSA32.txt')
calibrator = omni.RedundantCalibrator(64)

calibrator.antennaLocationTolerance = .1
calibrator.badAntenna = []
calibrator.badUBL = [0, 1, 2]

antlist=[[ -1.82119623e-02,  -1.17231687e-02 , -1.65467362e-02],
 [  3.98503466e+00,   9.63369375e-03 ,  2.37189551e-02],
 [  7.98828128e+00,   3.09905563e-02 ,  6.39846464e-02],
 [  1.19915279e+01,   5.23474188e-02,   1.04250338e-01],
 [ -5.14142894e-01,   1.20011739e+02,   3.95129958e-02],
 [  3.48910373e+00,   1.20033096e+02,   7.97786871e-02],
 [  7.49235035e+00,   1.20054452e+02,   1.20044378e-01],
 [  1.14955970e+01,   1.20075809e+02,   1.60310070e-01],
 [ -2.66177428e-01,   6.00000077e+01,   1.14831298e-02],
 [  3.73706919e+00,   6.00213646e+01,   5.17488211e-02],
 [  7.74031582e+00,   6.00427215e+01,   9.20145124e-02],
 [  1.17435624e+01,   6.00640783e+01,   1.32280204e-01],
 [ -7.62108360e-01,   1.80023470e+02,   6.75428618e-02],
 [  3.24113826e+00,   1.80044826e+02,   1.07808553e-01],
 [  7.24438488e+00,   1.80066183e+02,   1.48074244e-01],
 [  1.12476315e+01,   1.80087540e+02,   1.88339936e-01],
 [ -1.42194695e-01,   2.99941423e+01,  -2.53180318e-03],
 [  3.86105193e+00,   3.00154991e+01,   3.77338881e-02],
 [  7.86429855e+00,   3.00368560e+01,   7.79995794e-02],
 [  1.18675452e+01,   3.00582129e+01,   1.18265271e-01],
 [ -6.38125627e-01,   1.50017604e+02,   5.35279288e-02],
 [  3.36512099e+00,   1.50038961e+02,   9.37936201e-02],
 [  7.36836762e+00,   1.50060318e+02,   1.34059311e-01],
 [  1.13716142e+01,   1.50081675e+02,   1.74325003e-01],
 [ -3.90160161e-01,   9.00058732e+01,   2.54980628e-02],
 [  3.61308646e+00,   9.00272301e+01,   6.57637541e-02],
 [  7.61633308e+00,   9.00485869e+01,   1.06029445e-01],
 [  1.16195797e+01,   9.00699438e+01,   1.46295137e-01],
 [ -8.86091092e-01,   2.10029335e+02,   8.15577948e-02],
 [  3.11715553e+00,   2.10050692e+02,   1.21823486e-01],
 [  7.12040215e+00,   2.10072049e+02,   1.62089177e-01],
 [  1.11236488e+01,   2.10093406e+02,   2.02354869e-01]]
flat=[ele for bl in antlist for ele in bl]
calibrator.antennaLocation =np.reshape(np.array(flat),(len(flat)/3,3))
antloc=calibrator.antennaLocation

#nAntenna
nAntenna=len(calibrator.antennaLocation)


#find out UBL
##########################################################################################
#antloc has the form of a nested list with dimension nant*3, returns a np array of unique baselines
def UBL(antlist,precision):
	flat=[ele for bl in antlist for ele in bl]
	calibrator.antennaLocation =np.reshape(np.array(flat),(len(flat)/3,3))
	antloc=calibrator.antennaLocation
	ubllist=np.array([np.array([0,0,0])]);
	for i in range(len(antloc)):
		for j in range(i,len(antloc)):
			bool = False;
			for bl in ubllist:
				bool = bool or (np.linalg.norm(antloc[i]-antloc[j]-bl)<precision or np.linalg.norm(antloc[i]-antloc[j]+bl)<precision)
			if bool == False:			
				ubllist = np.concatenate((ubllist,[antloc[i]-antloc[j]]))
	ubllist=np.delete(ubllist,0,0)
	return ubllist


ubllist=np.zeros((1,3));
precision=calibrator.antennaLocationTolerance;
for i in range(len(antloc)):
	for j in range(i,len(antloc)):
		bool = False;
		for bl in ubllist:
			bool = bool or (np.linalg.norm(antloc[i]-antloc[j]-bl)<precision or np.linalg.norm(antloc[i]-antloc[j]+bl)<precision)
		if bool == False:			
			ubllist = np.concatenate((ubllist,[antloc[i]-antloc[j]]))
ubllist=np.delete(ubllist,0,0)
nUBL=len(ubllist)-len(calibrator.badUBL);

#################################################################################################
def dis(a1,a2):
	return np.linalg.norm(np.array(a1)-np.array(a2))


#find nBaseline (include auto baselines)
badbl=[ubllist[i] for i in calibrator.badUBL]
nbl=0;
for i in range(len(antloc)):
	for j in range(i,len(antloc)):
		bool=False
		for bl in badbl:
			bool=bool or dis(antloc[i]-antloc[j],bl)<precision or dis(antloc[i]-antloc[j],-bl)<precision
		if bool==False:
			nbl+=1
print nbl
		
