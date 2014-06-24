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

correctinfo = omni.read_redundantinfo('/home/ericy/omnical/redundant_info/redundantinfo_badant2_4_5_6_8_9_10badubl1_2_3_5_8_9.txt')
calibrator = omni.RedundantCalibrator(32)

calibrator.antennaLocationTolerance = .1
calibrator.badAntenna = [1,3,4,5,7,8,9]
calibrator.badUBL = [0, 1, 2,4,7,8]

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
calibrator.antennaLocation=np.reshape(np.array(flat),(len(flat)/3,3))

calibrator.compute_redundantinfo()

info1=correctinfo
info2=calibrator.info

#input two different redundant info, output True if they are the same and False if they are different
def compare_info(info1,info2,printfirstdiff=True):
	try:
		floatkeys=['antloc','ubl']
		intkeys = ['nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
		infomatrices=['A','B','At','Bt']
		allkeys=['antloc','ubl','nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB','A','B','At','Bt']
		diff=[]
		#10**5 for floating point errors
		for key in floatkeys:	
			diff.append(round(10**5*la.norm(info1[key]-info2[key]))==0)
		for key in intkeys:	
			diff.append(la.norm(info1[key]-info2[key])==0)
		for key in infomatrices:
			diff.append(la.norm((info1[key]-info2[key]).todense())==0)
		for i in info1['ublindex']-info2['ublindex']:
			diff.append(la.norm(i)==0)
		bool = True
		for i in diff:
			bool = bool and i
		#print the first key found different
		if printfirstdiff and bool == False:
			for i in range(len(diff)):
				if diff[i] == False:
					print allkeys[i]
					break
		return bool
	except ValueError:
		print "info doesn't have the same shape"
		return False
	
print compare_info(info1,info2)

print calibrator.get_baseline([1,2])


def compute_UBL(tolerance = 0.1):
	#check if the tolerance is not a string
	if type(tolerance) == str:
		raise Exception("tolerance needs to be number not string")
		return
	#remove the bad antennas
	nant=len(calibrator.antennaLocation)
	subsetant=[i for i in range(nant) if i not in calibrator.badAntenna]
	nAntenna=len(subsetant)
	antloc = np.array([calibrator.antennaLocation[ant] for ant in subsetant])
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

test=compute_UBL()
test2=calibrator.compute_UBL()



#old compute_UBL
	#def compute_UBL(self,tolerance = 0.1):
		##check if the tolerance is not a string
		#if type(tolerance) == str:
			#raise Exception("tolerance needs to be number not string")
			#return
		##remove the bad antennas
		#nant=len(self.antennaLocation)
		#subsetant=[i for i in range(nant) if i not in self.badAntenna]
		#nAntenna=len(subsetant)
		#antloc=[self.antennaLocation[ant] for ant in subsetant]
		#ubllist=np.array([np.array([0,0,0])]);
		#for i in range(len(antloc)):
			#for j in range(i+1,len(antloc)):
				#bool = False;
				#for bl in ubllist:
					#bool = bool or (la.norm(antloc[i]-antloc[j]-bl)<tolerance or la.norm(antloc[i]-antloc[j]+bl)<tolerance)
				#if bool == False:
					#ubllist = np.concatenate((ubllist,[antloc[j]-antloc[i]]))
		#ublall = np.delete(ubllist,0,0)
		#return ublall



'''
#nAntenna and subsetant : get rid of the bad antennas
nant=len(calibrator.antennaLocation)
subsetant=[i for i in range(nant) if i not in calibrator.badAntenna]
nAntenna=len(subsetant)
antloc=[calibrator.antennaLocation[ant] for ant in subsetant]
calibrator.totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(nant)])

#find out UBL
##########################################################################################
#antloc has the form of a nested list with dimension nant*3, returns a np array of unique baselines
def UBL(antloc,tolerance):
	ubllist=np.array([np.array([0,0,0])]);
	for i in range(len(antloc)):
		for j in range(i+1,len(antloc)):
			bool = False;
			for bl in ubllist:
				bool = bool or (np.linalg.norm(antloc[i]-antloc[j]-bl)<tolerance or np.linalg.norm(antloc[i]-antloc[j]+bl)<tolerance)
			if bool == False:			
				ubllist = np.concatenate((ubllist,[antloc[j]-antloc[i]]))
	ublall=np.delete(ubllist,0,0)
	return ublall


#ubllist=np.zeros((1,3));
#tolerance=calibrator.antennaLocationTolerance;
#for i in range(len(antloc)):
	#for j in range(i+1,len(antloc)):
		#bool = False;
		#for bl in ubllist:
			#bool = bool or (np.linalg.norm(antloc[i]-antloc[j]-bl)<tolerance or np.linalg.norm(antloc[i]-antloc[j]+bl)<tolerance)
		#if bool == False:			
			#ubllist = np.concatenate((ubllist,[antloc[j]-antloc[i]]))
			
tolerance=calibrator.antennaLocationTolerance;		
ublall=UBL(antloc,tolerance)
ubl=np.delete(ublall,np.array(calibrator.badUBL),0)
nUBL=len(ubl);

#################################################################################################
#calculate the norm of the difference of two vectors
def dis(a1,a2):
	return np.linalg.norm(np.array(a1)-np.array(a2))


#find nBaseline (include auto baselines) and subsetbl
badbl=[ublall[i] for i in calibrator.badUBL]
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

#from antenna pair to baseline index (pair[i,j], i>j)			
#def toBaseline(pair,n=nAntenna):
	#j=pair[0];
	#i=pair[1];	
	#return j*(j+1)/2+i
def toBaseline(pair,tvid=calibrator.totalVisibilityId):
	sortp=np.array(sorted(pair))
	for i in range(len(tvid)):
		if tvid[i][0] == sortp[0] and tvid[i][1] == sortp[1]:
			return i
	return 'no match'
subsetbl=np.array([toBaseline(bl,calibrator.totalVisibilityId) for bl in goodpairs])

##########################################
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
		
###########################################
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
		
###############################################
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

##################################################
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


################################################################
#bl1dmatrix: a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
		#I suppose 99999 for bad and auto baselines?
bl1dmatrix=99999*np.ones([nAntenna,nAntenna],dtype='int16')
for i in range(len(crosspair)):
	bl1dmatrix[crosspair[i][1]][crosspair[i][0]]=i
	bl1dmatrix[crosspair[i][0]][crosspair[i][1]]=i

################################################################
#degenM:
a=[] 
for i in range(len(antloc)):
	a.append(np.append(antloc[i],1))
a=np.array(a)
	
d=[]
for i in range(len(ubl)):
	d.append(np.append(ubl[i],0))
d=np.array(d)
	
	
m1=-a.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
m2=d.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
degenM = np.append(m1,m2,axis=0)

################################################################
#A: A matrix for logcal amplitude

A=np.zeros([len(crosspair),nAntenna+len(ubl)])
for i in range(len(crosspair)):
	A[i][crosspair[i][0]]=1
	A[i][crosspair[i][1]]=1
	A[i][nAntenna+bltoubl[i]]=1
A=sps.csr_matrix(A)

################################################################
#B: B matrix for logcal phase

B=np.zeros([len(crosspair),nAntenna+len(ubl)])
for i in range(len(crosspair)):
	B[i][crosspair[i][0]]=reverse[i]*1
	B[i][crosspair[i][1]]=reverse[i]*-1
	B[i][nAntenna+bltoubl[i]]=1
B=sps.csr_matrix(B)

#########################################
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
info['ncross']=ncross
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
		info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense())#(AtA)^-1
		info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense())#(BtB)^-1
		info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
		info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
		info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
		info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
		info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
		info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB

'''


