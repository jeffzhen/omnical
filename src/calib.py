'''XXX DOCSTRING'''
# XXX lots of imports... are all necessary?  can code be separated into files with smaller dependency lists?
# XXX this file has gotten huge. need to break into smaller files
# XXX clean house on commented code?
# XXX obey python style conventions
import math, random, traceback, ephem, string, commands, shutil, resource, threading, time
import multiprocessing as mp
from time import ctime
import aipy as ap
import struct
import numpy as np
import os, sys
import _omnical as _O
from info import RedundantInfo
import warnings
from array import array
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy as sp
    import scipy.sparse as sps
    import scipy.linalg as la
    import scipy.signal as ss
    import scipy.ndimage.filters as sfil
    from scipy import interpolate
    try:
        from numpy import nanmedian as nanmedian
    except:
        print "WARNING: using scipy's nanmedian function with is much slower than numpy.nanmedian. Consider numpy 1.9+."
        from scipy.stats import nanmedian

__version__ = '4.0.4'

julDelta = 2415020.# =julian date - pyephem's Observer date
PI = np.pi
TPI = 2 * np.pi

def apply_calpar(data, calpar, visibilityID):
    '''apply complex calpar for all antennas onto all baselines, calpar's dimension will be assumed to mean: 1D: constant over time and freq; 2D: constant over time; 3D: change over time and freq'''
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

def apply_calpar2(data, calpar, calpar2, visibilityID):
    '''apply complex calpar for all antennas onto all baselines, calpar's dimension will be assumed to mean: 1D: constant over time and freq; 2D: constant over time; 3D: change over time and freq'''
    METHODNAME = "*apply_calpar2*"
    if calpar.shape[-1] <= np.amax(visibilityID) or data.shape[-1] != len(visibilityID) or calpar.shape != calpar2.shape:
        raise Exception("Dimension mismatch! Either number of antennas in calpar " + str(calpar.shape[-1]) + " is less than implied in visibility ID "  + str(1 + np.amax(visibilityID)) + ", or the length of the last axis of data "  + str(data.shape[-1]) + " is not equal to length of visibilityID "  + str(len(visibilityID)) + ", or calpars have different dimensions:" + str(calpar.shape) + str(calpar.shape) + '.')
    if len(calpar.shape) == 3 and len(data.shape) == 3 and calpar.shape[:2] == data.shape[:2]:
        return data/(np.conjugate(calpar[:,:,visibilityID[:,0].astype(int)]) * calpar2[:,:,visibilityID[:,1].astype(int)])
    elif len(calpar.shape) == 2 and (len(data.shape) == 3 or len(data.shape) == 2) and calpar.shape[0] == data.shape[-2]:
        return data/(np.conjugate(calpar[:,visibilityID[:,0].astype(int)]) * calpar2[:,visibilityID[:,1].astype(int)])
    elif len(calpar.shape) == 1 and len(data.shape) <= 3:
        return data/(np.conjugate(calpar[visibilityID[:,0].astype(int)]) * calpar2[visibilityID[:,1].astype(int)])
    else:
        raise Exception("Dimension mismatch! I don't know how to interpret data dimension of " + str(data.shape) + " and calpar dimension of " + str(calpar.shape) + ".")

# XXX utility function, should be separate file
def stdmatrix(length, polydegree):
    '''to find out the error in fitting y by a polynomial poly(x), one compute error vector by (I-A.(At.A)^-1 At).y, where Aij = i^j. This function returns (I-A.(At.A)^-1 At)'''
    A = np.array([[i**j for j in range(polydegree + 1)] for i in range(length)], dtype='int')
    At = A.transpose()
    return np.identity(length) - A.dot(la.pinv(At.dot(A), cond = 10**(-6)).dot(At))

def omnical2omnigain(omnicalPath, utctimePath, info, outputPath = None):
    '''outputPath should be a path without extensions like .omnigain which will be appended'''
    if outputPath is None:
        outputPath = omnicalPath.replace('.omnical', '')

    #info = redundantCalibrator.info

    if not os.path.isfile(utctimePath):
        raise Exception("File %s does not exist!"%utctimePath)
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
    #opchisq[:, 3::2] = calpars[:, :, 0]#number of lincal iters
    opchisq[:, 3:] = calpars[:, :, 2]#chisq which is sum of squares of errors in each visbility

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

def _redcal(data, rawCalpar, Info, additivein, additive_out, removedegen=0, uselogcal=1, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit = 1, trust_period = 1):
    '''same as _O.redcal, but does not return additiveout. Rather it puts additiveout into an inputted container'''

    np_rawCalpar = np.frombuffer(rawCalpar, dtype='float32')
    np_rawCalpar.shape=(data.shape[0], data.shape[1], len(rawCalpar) / data.shape[0] / data.shape[1])
    #print np_rawCalpar.dtype, np_rawCalpar.shape

    np_additive_out = np.frombuffer(additive_out, dtype='complex64')
    np_additive_out.shape = data.shape
    _O.redcal2(data, np_rawCalpar, Info, additivein, np_additive_out, removedegen=removedegen, uselogcal=uselogcal, maxiter=int(maxiter), conv=float(conv), stepsize=float(stepsize), computeUBLFit = int(computeUBLFit), trust_period = int(trust_period))

    #np_additive_out = _O.redcal(data, np_rawCalpar, Info, additivein, removedegen=removedegen, uselogcal=uselogcal, maxiter=int(maxiter), conv=float(conv), stepsize=float(stepsize), computeUBLFit = int(computeUBLFit), trust_period = int(trust_period))
    #additive_out[:len(additive_out)/2] = np.real(np_additive_out.flatten())
    #additive_out[len(additive_out)/2:] = np.imag(np_additive_out.flatten())


#  ___        _              _          _    ___      _ _ _             _           
# | _ \___ __| |_  _ _ _  __| |__ _ _ _| |_ / __|__ _| (_) |__ _ _ __ _| |_ ___ _ _ 
# |   / -_) _` | || | ' \/ _` / _` | ' \  _| (__/ _` | | | '_ \ '_/ _` |  _/ _ \ '_|
# |_|_\___\__,_|\_,_|_||_\__,_\__,_|_||_\__|\___\__,_|_|_|_.__/_| \__,_|\__\___/_|  

class RedundantCalibrator:
    '''This class is the main tool for performing redundant calibration on data sets. 
    For a given redundant configuration, say 32 antennas with 3 bad antennas, the 
    user should create one instance of Redundant calibrator and reuse it for all data 
    collected from that array. In general, upon creating an instance, the user need 
    to create the info field of the instance by either computing it or reading it 
    from a text file. readyForCpp(verbose = True) should be a very helpful function 
    to provide information on what information is missing for running the calibration.'''
    def __init__(self, nTotalAnt, info = None):
        methodName = '.__init__.'
        self.className = '.RedundantCalibrator.'
        self.nTotalAnt = nTotalAnt
        self.nTotalBaselineAuto = (self.nTotalAnt + 1) * self.nTotalAnt / 2
        self.nTotalBaselineCross = (self.nTotalAnt - 1) * self.nTotalAnt / 2
        self.antennaLocation = np.zeros((self.nTotalAnt, 3))
        side = int(nTotalAnt**.5)
        for a in range(nTotalAnt):
            self.antennaLocation[a] = np.array([a/side, a%side, 0])
        self.antennaLocationTolerance = 10**(-6)
        self.badAntenna = []
        self.badUBL = []
        self.badUBLpair = []
        self.ubl2goodubl = None
        self.totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(self.nTotalAnt)])#PAPER miriad convention by default
        self.totalVisibilityId_dic = None
        self.totalVisibilityUBL = None
        self.Info = None
        self.removeDegeneracy = True
        self.removeAdditive = False
        self.removeAdditivePeriod = -1
        self.convergePercent = 0.01 #convergence criterion in relative change of chi^2. By default it stops when reaches 0.01, namely 1% decrease in chi^2.
        self.maxIteration = 50 #max number of iterations in lincal
        self.stepSize = 0.3 #step size for lincal. (0, 1]. < 0.4 recommended.
        self.computeUBLFit = True
        self.trust_period = 1 #How many time slices does lincal start from logcal result rather than the previous time slice's lincal result. default 1 means always start from logcal. if 10, it means lincal start from logcal results (or g = 1's) every 10 time slices

        self.nTime = 0
        self.nFrequency = 0

        self.utctime = None
        self.rawCalpar = None
        self.omnichisq = None
        self.omnigain = None
        self.omnifit = None

        if info is not None: # XXX what is the point of leaving info == None?
            if type(info) == str:
                self.read_redundantinfo(info)
            else:
                self.Info = info
    def __repr__(self,): return self.__str__()
    def __str__(self,):
        if self.Info is None:
            return "<Uninitialized %i antenna RedundantCalibrator with no RedundantInfo.>"%self.nTotalAnt
        else:
            return "<RedundantCalibrator for an %i antenna array: %i good baselines including %i good antennas and %i unique baselines.>"%(self.nTotalAnt, len(self.Info.crossindex), self.Info.nAntenna, self.Info.nUBL)

    def __getattr__(self, name): # XXX why not just inherit from info if accessing all attributes of info?
        try: return self.Info.__getattribute__(name)
        except: raise AttributeError("RedundantCalibrator has no attribute named %s"%name)

    def read_redundantinfo(self, filename, txtmode=False, verbose=False):
        '''redundantinfo is necessary for running redundant calibration. The text file 
        should contain 29 lines each describes one item in the info.'''
        self.Info = RedundantInfo(filename=filename, txtmode=txtmode, verbose=verbose)
        self.totalVisibilityId = self.Info.totalVisibilityId # XXX might this raise an exception?
        #try: self.totalVisibilityId = self.Info.totalVisibilityId
        #except(KeyError): self.Info.totalVisibilityId = self.totalVisibilityId

    def write_redundantinfo(self, filename, overwrite=False, verbose=False):
        self.Info.tofile(filename, overwrite=overwrite, verbose=verbose)

    def read_arrayinfo(self, arrayinfopath, verbose = False):
        '''array info is the minimum set of information to uniquely describe a 
        redundant array, and is needed to compute redundant info. It includes, 
        in each line, bad antenna indices, bad unique baseline indices, tolerance 
        of error when checking redundancy, antenna locations, and visibility's 
        antenna pairing conventions. Unlike redundant info which is a self-contained 
        dictionary, items in array info each have their own fields in the instance.'''
        methodName = ".read_arrayinfo."
        if not os.path.isfile(arrayinfopath):
            raise IOError(self.className + methodName + "Error: Array info file " + arrayinfopath + " doesn't exist!")
        with open(arrayinfopath) as f:
            rawinfo = [[float(x) for x in line.split()] for line in f]
        if verbose:
            print self.className + methodName + " MSG:",  "Reading", arrayinfopath, "...",

        self.badAntenna = np.array(rawinfo[0]).astype(int)
        if self.badAntenna[0] < 0:
            self.badAntenna = np.zeros(0)

        if len(np.array(rawinfo[1])) == 0 or len(np.array(rawinfo[1]))%2 != 0 or min(np.array(rawinfo[1]).astype(int)) < 0:
            self.badUBLpair = np.array([])
            #raise Exception(self.className + methodName +"Error: Format error in " + arrayinfopath + "badUBL should be specified by pairs of antenna, not odd numbers of antenna")
        else:
            rawpair = np.array(rawinfo[1]).astype(int)
            self.badUBLpair = np.reshape(rawpair,(len(rawpair)/2,2))
        #if self.badUBL[0] < 0:#todonow
        #    self.badUBL = np.zeros(0)

        self.antennaLocationTolerance = rawinfo[2][0]

        for a in range(len(self.antennaLocation)):
            if len(rawinfo[a + 3]) != 3:
                raise ValueError(self.className + methodName + "Error: Format error in " + arrayinfopath + ": The antenna locations should start on the 4th line, with 3 numbers in each line!")
            else:
                self.antennaLocation[a] = np.array(rawinfo[a + 3])

        #for bl in range(len(self.totalVisibilityId)):
            #if len(rawinfo[bl + 3 + len(self.antennaLocation)]) != 2:
                #raise Exception(self.className + methodName + "Error: Format error in " + arrayinfopath + ": The baseline to antenna mapping should start after antenna locations, with 2 numbers (conj index, index) in each line!")
            #else:
        bl = 0
        self.totalVisibilityId = []
        max_bl_cnt = self.nTotalAnt * (self.nTotalAnt + 1) / 2
        maxline = len(rawinfo)
        while len(rawinfo[bl + 3 + len(self.antennaLocation)]) == 2:
            if bl >= max_bl_cnt:
                raise Exception("Number of total visibility ids exceeds the maximum possible number of baselines of %i"%(max_bl_cnt))
            self.totalVisibilityId.append(np.array(rawinfo[bl + 3 + len(self.antennaLocation)]).astype(int))
            bl = bl + 1
            if bl + 3 + len(self.antennaLocation) >= maxline:
                break
        self.totalVisibilityId = np.array(self.totalVisibilityId).astype(int)
        if verbose:
            print "Total number of visibilities:", bl,
            print "Bad antenna indices:", self.badAntenna,
            print "Bad UBL indices:", self.badUBLpair


    def lincal(self, data, additivein, nthread = None, verbose = False):
        '''for best performance, try setting nthread to larger than number of cores.'''
        if data.ndim != 3 or data.shape[-1] != len(self.totalVisibilityId):
            raise ValueError("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline(%i)"%len(self.totalVisibilityId))
        if data.shape != additivein.shape:
            raise ValueError("Data shape error: data and additive in have different shapes.")
        self.nTime = len(data)
        self.nFrequency = len(data[0])
        if self.rawCalpar is None:
            self.rawCalpar = np.zeros((self.nTime, self.nFrequency, 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)), dtype = 'float32')
        elif self.rawCalpar.shape != (len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)):
            raise ValueError("ERROR: lincal called without a properly shaped self.rawCalpar! Excpeted shape is (%i, %i, %i)!"%(len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)))
        if nthread is None:
            nthread = nthread = min(mp.cpu_count() - 1, self.nFrequency)
        if nthread < 2:
            return _O.redcal(data[:,:,self.Info.subsetbl], self.rawCalpar, self.Info, additivein[:,:,self.Info.subsetbl], removedegen = int(self.removeDegeneracy), uselogcal = 0, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit), trust_period = self.trust_period)
        else:
            return self._redcal_multithread(data, additivein, 0, nthread, verbose = verbose)        ##self.chisq = self.rawCalpar[:, :, 2]
        ##self.calpar = np.zeros((len(self.rawCalpar), len(self.rawCalpar[0]), self.nTotalAnt), dtype='complex64')
        ##self.calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
        ##self.bestfit = self.rawCalpar[:, :, (3 + 2 * self.Info.nAntenna):: 2] + 1.j * self.rawCalpar[:, :, (4 + 2 * self.Info.nAntenna):: 2]

    def logcal(self, data, additivein, nthread = None, verbose = False):
        '''XXX DOCSTRING'''
        if data.ndim != 3 or data.shape[-1] != len(self.totalVisibilityId):
            raise ValueError("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline(%i)"%len(self.totalVisibilityId))
        if data.shape != additivein.shape:
            raise ValueError("Data shape error: data and additive in have different shapes.")
        self.nTime = len(data)
        self.nFrequency = len(data[0])
        self.rawCalpar = np.zeros((len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)), dtype = 'float32')

        if nthread is None:
            nthread = min(mp.cpu_count() - 1, self.nFrequency)
        if nthread < 2:
            return _O.redcal(data[:,:,self.Info.subsetbl], self.rawCalpar, self.Info, additivein[:,:,self.Info.subsetbl], removedegen = int(self.removeDegeneracy), uselogcal = 1, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit))
        else:
            return self._redcal_multithread(data, additivein, 1, nthread, verbose = verbose)

    def _redcal_multithread(self, data, additivein, uselogcal, nthread, verbose = False):
        '''XXX DOCSTRING'''
        #if data.ndim != 3 or data.shape[-1] != len(self.totalVisibilityId):
            #raise ValueError("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline(%i)"%len(self.totalVisibilityId))
        #if data.shape != additivein.shape:
            #raise ValueError("Data shape error: data and additive in have different shapes.")
        #self.nTime = len(data)
        #self.nFrequency = len(data[0])
        #self.rawCalpar = np.zeros((len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)), dtype = 'float32')
        nthread = min(nthread, self.nFrequency)
        additiveouts = {}
        np_additiveouts = {}
        #additiveout = np.zeros_like(data[:, :, self.Info.subsetbl])
        rawCalpar = {}
        np_rawCalpar = {}
        threads = {}
        fchunk = {}
        chunk = int(self.nFrequency) / int(nthread)
        excess = int(self.nFrequency) % int(nthread)
        kwarg = {"removedegen": int(self.removeDegeneracy), "uselogcal": uselogcal, "maxiter": int(self.maxIteration), "conv": float(self.convergePercent), "stepsize": float(self.stepSize), "computeUBLFit": int(self.computeUBLFit), "trust_period": self.trust_period}

        for i in range(nthread):
            if excess == 0:
                fchunk[i] = (i * chunk, min((1 + i) * chunk, self.nFrequency),)
            elif i < excess:
                fchunk[i] = (i * (chunk+1), min((1 + i) * (chunk+1), self.nFrequency),)
            else:
                fchunk[i] = (fchunk[i-1][1], min(fchunk[i-1][1] + chunk, self.nFrequency),)
            #if verbose:
                #print fchunk[i],
            rawCalpar[i] = mp.RawArray('f', self.nTime * (fchunk[i][1] - fchunk[i][0]) * (self.rawCalpar.shape[2]))
            np_rawCalpar[i] = np.frombuffer(rawCalpar[i], dtype='float32')
            np_rawCalpar[i].shape = (self.rawCalpar.shape[0], fchunk[i][1]-fchunk[i][0], self.rawCalpar.shape[2])
            np_rawCalpar[i][:] = self.rawCalpar[:, fchunk[i][0]:fchunk[i][1]]

            additiveouts[i] = mp.RawArray('f', self.nTime * (fchunk[i][1] - fchunk[i][0]) * len(self.Info.subsetbl) * 2)#factor of 2 for re/im
            np_additiveouts[i] = np.frombuffer(additiveouts[i], dtype='complex64')
            np_additiveouts[i].shape = (data.shape[0], fchunk[i][1]-fchunk[i][0], len(self.Info.subsetbl))

            threads[i] = mp.Process(target = _redcal, args = (data[:, fchunk[i][0]:fchunk[i][1], self.Info.subsetbl], rawCalpar[i], self.Info, additivein[:, fchunk[i][0]:fchunk[i][1], self.Info.subsetbl], additiveouts[i]), kwargs=kwarg)
            #threads[i] = threading.Thread(target = _O.redcal2, args = (data[:, fchunk[i][0]:fchunk[i][1], self.Info.subsetbl], rawCalpar[i], self.Info, additivein[:, fchunk[i][0]:fchunk[i][1], self.Info.subsetbl], additiveouts[i]), kwargs=kwarg)
        if verbose:
            print "Starting %s Process"%cal_name[uselogcal],
            sys.stdout.flush()
        for i in range(nthread):
            if verbose:
                print "#%i"%i,
                sys.stdout.flush()
            threads[i].start()
        if verbose:
            print "Finished Process",
        for i in range(nthread):
            threads[i].join()
            if verbose:
                print "#%i"%i,
        if verbose:
            print ""
            sys.stdout.flush()
        self.rawCalpar = np.concatenate([np_rawCalpar[i] for i in range(nthread)],axis=1)
        return np.concatenate([np_additiveouts[i] for i in range(nthread)],axis=1)

    def get_calibrated_data(self, data, additivein = None):
        '''XXX DOCSTRING'''
        if data.ndim != 3 or data.shape != (self.nTime, self.nFrequency, len(self.totalVisibilityId)):
            raise ValueError("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline (%i, %i, %i)"%(self.nTime, self.nFrequency, len(self.totalVisibilityId)))
        if additivein is not None and data.shape != additivein.shape:
            raise ValueError("Data shape error: data and additivein have different shapes.")
        if data.shape[:2] != self.rawCalpar.shape[:2]:
            raise ValueError("Data shape error: data and self.rawCalpar have different first two dimensions.")

        calpar = np.ones((len(self.rawCalpar), len(self.rawCalpar[0]), self.nTotalAnt), dtype='complex64')
        calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
        if additivein is None:
            return apply_calpar(data, calpar, self.totalVisibilityId)
        else:
            return apply_calpar(data - additivein, calpar, self.totalVisibilityId)

    def get_modeled_data(self):
        '''XXX DOCSTRING'''
        if self.rawCalpar is None:
            raise ValueError("self.rawCalpar doesn't exist. Please calibrate first using logcal() or lincal().")
        if len(self.totalVisibilityId) <= np.max(self.Info.subsetbl):
            raise ValueError("self.totalVisibilityId of length %i is shorter than max index in subsetbl %i. Probably you are using an outdated version of redundantinfo."%(len(self.totalVisibilityId), np.max(self.Info.subsetbl)))
        mdata = np.zeros((self.rawCalpar.shape[0], self.rawCalpar.shape[1], len(self.totalVisibilityId)), dtype='complex64')
        mdata[..., self.Info.subsetbl[self.Info.crossindex]] = (self.rawCalpar[..., 3 + 2 * (self.Info.nAntenna)::2] + 1.j * self.rawCalpar[..., 4 + 2 * (self.Info.nAntenna)::2])[..., self.Info.bltoubl]
        mdata[..., self.Info.subsetbl[self.Info.crossindex]] = np.abs(mdata[..., self.Info.subsetbl[self.Info.crossindex]]) * np.exp(self.Info.reversed * 1.j * np.angle(mdata[..., self.Info.subsetbl[self.Info.crossindex]])) * 10.**(self.rawCalpar[..., 3 + self.Info.bl2d[self.Info.crossindex,0]] + self.rawCalpar[..., 3 + self.Info.bl2d[self.Info.crossindex,1]]) * np.exp(-1.j * self.rawCalpar[..., 3 + self.Info.nAntenna + self.Info.bl2d[self.Info.crossindex,0]] + 1.j * self.rawCalpar[..., 3 + self.Info.nAntenna + self.Info.bl2d[self.Info.crossindex,1]])
        return mdata


    def get_omnichisq(self):
        '''XXX DOCSTRING'''
        if self.utctimes is None or self.rawCalpar is None:
            raise Exception("Error: either self.utctimes or self.rawCalpar does not exist.")
        if len(self.utctimes) != len(self.rawCalpar):
            raise Exception("Error: length of self.utctimes is not equal to self.rawCalpar. One of them is wrong.")
        jd = np.zeros((len(self.utctimes), 2), dtype='float32')#Julian dat is the only double in this whole thing so im storing it in two chunks as float
        sa = ephem.Observer()
        for utctime, t in zip(self.utctimes, range(len(self.utctimes))):
            sa.date = utctime
            jd[t, :] = struct.unpack('ff', struct.pack('d', sa.date + julDelta))

        omnichisq = np.zeros((self.nTime, 2 + 1 + self.nFrequency), dtype = 'float32')
        omnichisq[:, :2] = jd
        omnichisq[:, 2] = float(self.nFrequency)
        omnichisq[:, 3:] = self.rawCalpar[:, :, 2]#chisq which is sum of squares of errors in each visbility
        return omnichisq

    def get_omnigain(self):
        '''XXX DOCSTRING'''
        if self.utctimes is None or self.rawCalpar is None:
            raise Exception("Error: either self.utctimes or self.rawCalpar does not exist.")
        if len(self.utctimes) != len(self.rawCalpar):
            raise Exception("Error: length of self.utctimes is not equal to self.rawCalpar. One of them is wrong.")
        jd = np.zeros((len(self.utctimes), 2), dtype='float32')#Julian dat is the only double in this whole thing so im storing it in two chunks as float
        sa = ephem.Observer()
        for utctime, t in zip(self.utctimes, range(len(self.utctimes))):
            sa.date = utctime
            jd[t, :] = struct.unpack('ff', struct.pack('d', sa.date + julDelta))
        omnigain = np.zeros((self.nTime, self.Info.nAntenna, 2 + 1 + 1 + 2 * self.nFrequency), dtype = 'float32')
        omnigain[:, :, :2] = jd[:, None]
        omnigain[:, :, 2] = np.array(self.Info.subsetant).astype('float32')
        omnigain[:, :, 3] = float(self.nFrequency)
        gains = (10**self.rawCalpar[:, :, 3:(3 + self.Info.nAntenna)] * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna):(3 + 2 * self.Info.nAntenna)])).transpose((0,2,1))
        omnigain[:, :, 4::2] = np.real(gains)
        omnigain[:, :, 5::2] = np.imag(gains)
        return omnigain

    def get_omnifit(self):
        '''XXX DOCSTRING'''
        if self.utctimes is None or self.rawCalpar is None:
            raise Exception("Error: either self.utctimes or self.rawCalpar does not exist.")
        if len(self.utctimes) != len(self.rawCalpar):
            raise Exception("Error: length of self.utctimes is not equal to self.rawCalpar. One of them is wrong.")
        jd = np.zeros((len(self.utctimes), 2), dtype='float32')#Julian dat is the only double in this whole thing so im storing it in two chunks as float
        sa = ephem.Observer()
        for utctime, t in zip(self.utctimes, range(len(self.utctimes))):
            sa.date = utctime
            jd[t, :] = struct.unpack('ff', struct.pack('d', sa.date + julDelta))
        omnifit = np.zeros((self.nTime, self.Info.nUBL , 2 + 3 + 1 + 2 * self.nFrequency), dtype = 'float32')
        omnifit[:, :, :2] = jd[:, None]
        omnifit[:, :, 2:5] = np.array(self.Info.ubl).astype('float32')
        omnifit[:, :, 5] = float(self.nFrequency)
        omnifit[:, :, 6::2] = self.rawCalpar[:, :, 3 + 2 * self.Info.nAntenna::2].transpose((0,2,1))
        omnifit[:, :, 7::2] = self.rawCalpar[:, :, 3 + 2 * self.Info.nAntenna + 1::2].transpose((0,2,1))
        return omnifit

    def set_badUBL(self, badUBL):
        '''XXX DOCSTRING'''
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



    def diagnose(self, data = None, additiveout = None, flag = None, verbose = True, healthbar = 2, ubl_healthbar = 50, warn_low_redun = False, ouput_txt = False):
        '''XXX DOCSTRING'''
        errstate = np.geterr()
        np.seterr(invalid = 'ignore')

        if self.rawCalpar is None:
            raise Exception("No calibration has been performed since rawCalpar does not exist.")

        if flag is None:
            flag = np.zeros(self.rawCalpar.shape[:2], dtype='bool')
        elif flag.shape != self.rawCalpar.shape[:2]:
            raise TypeError('flag and self.rawCalpar have different shapes %s %s.'%(flag.shape, self.rawCalpar.shape[:2]))

        checks = 1
        timer = Timer()
        bad_count = np.zeros((3,self.Info.nAntenna), dtype='int')
        bad_ubl_count = np.zeros(self.Info.nUBL, dtype='int')
        median_level = nanmedian(nanmedian(self.rawCalpar[:,:,3:3+self.Info.nAntenna], axis= 0), axis= 1)
        bad_count[0] = np.array([(np.abs(self.rawCalpar[:,:,3+a] - median_level) >= .15)[~flag].sum() for a in range(self.Info.nAntenna)])**2
        #timer.tick(1)


        if data is not None and data.shape[:2] == self.rawCalpar.shape[:2]:
            checks += 1
            subsetbl = self.Info.subsetbl
            crossindex = self.Info.crossindex
            ncross = len(self.Info.crossindex)
            bl1dmatrix = self.Info.bl1dmatrix
            ant_level = np.array([np.median(np.abs(data[:, :, [subsetbl[crossindex[bl]] for bl in bl1dmatrix[a] if (bl < ncross and bl >= 0)]]), axis = -1) for a in range(self.Info.nAntenna)])
            #timer.tick(2)
            median_level = nanmedian(ant_level, axis = 0)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",category=RuntimeWarning)
                bad_count[1] = np.array([(np.abs(ant_level[a] - median_level)/median_level >= .667)[~flag].sum() for a in range(self.Info.nAntenna)])**2
        #timer.tick(2)

        if additiveout is not None and additiveout.shape[:2] == self.rawCalpar.shape[:2]:
            checks += 1

            subsetbl = self.Info.subsetbl
            crossindex = self.Info.crossindex
            ncross = len(self.Info.crossindex)
            bl1dmatrix = self.Info.bl1dmatrix
            ant_level = np.array([np.median(np.abs(additiveout[:, :, [crossindex[bl] for bl in bl1dmatrix[a] if bl < ncross]]), axis = 2) for a in range(self.Info.nAntenna)])
            #timer.tick(3)
            median_level = np.median(ant_level, axis = 0)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",category=RuntimeWarning)
                bad_count[2] = np.array([(np.abs(ant_level[a] - median_level)/median_level >= .667)[~flag].sum() for a in range(self.Info.nAntenna)])**2
            #timer.tick(3)
            ublindex = [np.array(index).astype('int')[:,2] for index in self.Info.ublindex]
            ubl_level = np.array([np.median(np.abs(additiveout[:, :, [crossindex[bl] for bl in ublindex[u]]]), axis = 2) for u in range(self.Info.nUBL)])
            median_level = np.median(ubl_level, axis = 0)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",category=RuntimeWarning)
                bad_ubl_count += np.array([((ubl_level[u] - median_level)/median_level >= .667)[~flag].sum() for u in range(self.Info.nUBL)])**2
            #print median_level
        #timer.tick(3)

        np.seterr(invalid = errstate['invalid'])
        bad_count = (np.mean(bad_count,axis=0)/float(np.sum(~flag))**2 * 100).astype('int')
        bad_ubl_count = (bad_ubl_count/float(self.nTime * self.nFrequency)**2 * 100).astype('int')
        if verbose:
            #print bad_ant_cnt, bad_ubl_cnt
            print "DETECTED BAD ANTENNA ABOVE HEALTH THRESHOLD %i: "%healthbar
            for a in range(len(bad_count)):
                if bad_count[a] > healthbar:
                    print "antenna #%i, vector = %s, badness = %i"%(self.Info.subsetant[a], self.Info.antloc[a], bad_count[a])
            #print ""
            if additiveout is not None and additiveout.shape[:2] == self.rawCalpar.shape[:2] and ubl_healthbar != 100:
                print "DETECTED BAD BASELINE TYPE ABOVE HEALTH THRESHOLD %i: "%ubl_healthbar
                for a in range(len(bad_ubl_count)):
                    if bad_ubl_count[a] > ubl_healthbar and (self.Info.ublcount[a] > 5 or (warn_low_redun)):
                        print "index #%i, vector = %s, redundancy = %i, badness = %i"%(a, self.Info.ubl[a], self.Info.ublcount[a], bad_ubl_count[a])
                #print ""
        if not ouput_txt:
            return bad_count, bad_ubl_count
        else:
            txt = ''
            txt += "DETECTED BAD ANTENNA ABOVE HEALTH THRESHOLD %i: \n"%healthbar
            for a in range(len(bad_count)):
                if bad_count[a] > healthbar:
                    txt += "antenna #%i, vector = %s, badness = %i\n"%(self.Info.subsetant[a], self.Info.antloc[a], bad_count[a])
            #print ""
            if additiveout is not None and additiveout.shape[:2] == self.rawCalpar.shape[:2] and ubl_healthbar != 100:
                txt += "DETECTED BAD BASELINE TYPE ABOVE HEALTH THRESHOLD %i: \n"%ubl_healthbar
                for a in range(len(bad_ubl_count)):
                    if bad_ubl_count[a] > ubl_healthbar and (self.Info.ublcount[a] > 5 or (warn_low_redun)):
                        txt += "index #%i, vector = %s, redundancy = %i, badness = %i\n"%(a, self.Info.ubl[a], self.Info.ublcount[a], bad_ubl_count[a])
            return txt

    def flag(self, mode = '12', twindow = 5, fwindow = 5, nsigma = 4, _dbg_plotter = None, _niter = 3):
        '''return true if flagged False if good and unflagged'''
        if self.rawCalpar is None or (self.rawCalpar[:,:,2] == 0).all():
            raise Exception("flag cannot be run before lincal.")

        chisq = np.copy(self.rawCalpar[:,:,2])
        nan_flag = np.isnan(np.sum(self.rawCalpar,axis=-1))|np.isinf(np.sum(self.rawCalpar,axis=-1))

        #chisq flag: spike_flag
        spike_flag = np.zeros_like(nan_flag)
        if '1' in mode:
            median_level = nanmedian(nanmedian(chisq))

            thresh = nsigma * (2. / (len(self.subsetbl) - self.nAntenna - self.nUBL + 2))**.5 # relative sigma is sqrt(2/k)

            for i in range(_niter):
                chisq[nan_flag|spike_flag] = 1e6 * median_level

                if twindow >= self.nTime:
                    filtered_tdir = np.ones(self.nTime)#will rescale anyways * np.min(np.median(chisq, axis = 1))
                else:
                    filtered_tdir = sfil.minimum_filter(np.median(chisq, axis = 1), size = twindow, mode='reflect')

                if fwindow >= self.nFrequency:
                    filtered_fdir = np.ones(self.nFrequency)#will rescale anyways * np.min(np.median(chisq, axis = 0))
                else:
                    filtered_fdir = sfil.minimum_filter(np.median(chisq, axis = 0), size = fwindow, mode='reflect')

                smoothed_chisq = np.outer(filtered_tdir, filtered_fdir)

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",category=RuntimeWarning)
                    smoothed_chisq = smoothed_chisq * np.median(chisq[~(nan_flag|spike_flag)] / smoothed_chisq[~(nan_flag|spike_flag)])

                    del_chisq = chisq - smoothed_chisq
                    del_chisq[(nan_flag|spike_flag)] = np.nan
                    estimate_chisq_sigma = np.nanstd(del_chisq,axis=0)
                    estimate_chisq_sigma[np.isnan(estimate_chisq_sigma)] = 0
                spike_flag = spike_flag | (np.abs(chisq - smoothed_chisq) >= np.minimum(smoothed_chisq * thresh, estimate_chisq_sigma * nsigma)) | (chisq == 0)

            if _dbg_plotter is not None:
                _dbg_plotter.imshow(np.abs(chisq - smoothed_chisq)/smoothed_chisq, vmin=-thresh, vmax=thresh, interpolation='none')
        #baseline fit flag
        if '2' in mode:
            nubl = 10
            short_ubl_index = np.argsort(np.linalg.norm(self.ubl, axis=1))[:min(nubl, self.nUBL)]
            shortest_ubl_vis = self.rawCalpar[:,:,3+2*self.nAntenna+2*short_ubl_index] + 1.j * self.rawCalpar[:,:,3+2*self.nAntenna+2*short_ubl_index+1]
            change_rate = np.median(np.abs(shortest_ubl_vis[:-1] - shortest_ubl_vis[1:]), axis = 2)
            nan_mask2 = np.isnan(change_rate)|np.isinf(change_rate)
            change_rate[nan_mask2] = 0

            if twindow >= self.nTime:
                filtered_tdir = np.ones(self.nTime - 1)#will rescale anyways * np.min(np.median(chisq, axis = 1))
            else:
                filtered_tdir = sfil.minimum_filter(np.median(change_rate, axis = 1), size = twindow, mode='reflect')

            if fwindow >= self.nFrequency:
                filtered_fdir = np.ones(self.nFrequency)#will rescale anyways * np.min(np.median(chisq, axis = 0))
            else:
                filtered_fdir = sfil.minimum_filter(np.median(change_rate, axis = 0), size = fwindow, mode='reflect')

            smoothed_change_rate = np.outer(filtered_tdir, filtered_fdir)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",category=RuntimeWarning)
                smoothed_change_rate = smoothed_change_rate * np.median(change_rate[~nan_mask2] / smoothed_change_rate[~nan_mask2])

            ubl_flag_short = (change_rate > 2 * smoothed_change_rate) | nan_mask2
            ubl_flag = np.zeros_like(spike_flag)
            ubl_flag[:-1] = ubl_flag_short
            ubl_flag[1:] = ubl_flag[1:] | ubl_flag_short
        else:
            ubl_flag = np.zeros_like(nan_flag)

        return_flag = (nan_flag|spike_flag|ubl_flag)
        return return_flag

    def compute_redundantinfo(self, arrayinfoPath = None, verbose = False, badAntenna = [], badUBLpair = [], antennaLocationTolerance = 1e-6):
        '''XXX DOCSTRING'''
        self.antennaLocationTolerance = antennaLocationTolerance
        self.badAntenna += badAntenna
        self.badUBLpair += badUBLpair

        if arrayinfoPath is not None and os.path.isfile(arrayinfoPath):
            self.read_arrayinfo(arrayinfoPath)
        if np.linalg.norm(self.antennaLocation) == 0:
            raise Exception("Error: compute_redundantinfo() called before self.antennaLocation is specified. Use configFilePath option when calling compute_redundantinfo() to specify array info file, or manually set self.antennaLocation for the RedundantCalibrator instance.")
        if verbose:
            timer = Timer()

        #antennalocation quality check: make sure there's no crazy constant added to antlocs
        bad_ant_mask = np.array([a in self.badAntenna for a in range(len(self.antennaLocation))]).astype('bool')
        array_center = la.norm(np.mean(self.antennaLocation[~bad_ant_mask], axis = 0))
        array_std = la.norm(np.std(self.antennaLocation[~bad_ant_mask], axis = 0))
        #print array_center, array_std
        if array_std / array_center < 1e-3:
            raise TypeError("Average antenna location is %s whereas the typical variation among locations is %s, which is too small and will cause many problems. Please remove the large overall offset from antenna locations."%(np.mean(self.antennaLocation[~bad_ant_mask], axis = 0), np.std(self.antennaLocation[~bad_ant_mask], axis = 0)))

        #nAntenna and subsetant : get rid of the bad antennas
        nant=len(self.antennaLocation)
        #subsetant=[i for i in range(nant) if i not in self.badAntenna]
        ant2goodant = -np.ones(len(self.antennaLocation), dtype=int)
        subsetant = []
        for a in range(len(self.antennaLocation)):
            if a not in self.badAntenna:
                subsetant.append(a)
                ant2goodant[a] = len(subsetant) - 1

        nAntenna=len(subsetant)
        antloc=[self.antennaLocation[ant] for ant in subsetant]
        if verbose:
            timer.tick('a')
        ##########################################################################################
        #find out ubl
        #use the function compute_UBL to find the ubl
        tolerance=self.antennaLocationTolerance;
        ublall=self.compute_UBL(tolerance)
        if verbose:
            timer.tick('b')
        #################################################################################################
        #calculate the norm of the difference of two vectors (just la.norm actually)
        def dis(a1,a2):
            return np.linalg.norm(np.array(a1)-np.array(a2))
        #find badUBL with badUBLpair
        def find_ublindex_all(pair):
            #print pair
            for i in range(len(ublall)):
                if dis(self.antennaLocation[pair[0]]-self.antennaLocation[pair[1]],ublall[i]) < tolerance or dis(self.antennaLocation[pair[0]]-self.antennaLocation[pair[1]],-ublall[i]) < tolerance:
                    return i
            return None
            #raise Exception("Error: something wrong in identifying badUBL from badUBLpair")    #delete this line afterwards
        #print self.badUBLpair, len(self.badUBLpair),self.badUBLpair[0]
        for p in self.badUBLpair:
            self.badUBL.append(find_ublindex_all(p))
        self.badUBL = [i for i in self.badUBL if i is not None]
        self.ubl2goodubl = -np.ones(len(ublall), dtype=int)
        goodu = 0
        for u in range(len(ublall)):
            if u not in self.badUBL:
                self.ubl2goodubl[u] = goodu
                goodu = goodu + 1

        #delete the bad ubl's
        ubl=np.delete(ublall,np.array(self.badUBL).astype('int'),0)
        nUBL=len(ubl);
        badbl=[ublall[i] for i in self.badUBL]
        #find nBaseline (include auto baselines) and subsetbl
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
        ####for a1, a2 in self.totalVisibilityId:
            ####i = ant2goodant[a1]
            ####j = ant2goodant[a2]
            ####if i >= 0 and j >= 0:
                ####bool=False
                ####for bl in badbl:
                    ####bool = bool or dis(antloc[i]-antloc[j],bl)<tolerance or dis(antloc[i]-antloc[j],-bl)<tolerance
                ####if bool == False:
                    ####nbl+=1
                    ####goodpairs.append([i,j])

        self.totalVisibilityId_dic = {}
        for bll, (a1, a2) in enumerate(self.totalVisibilityId):
            self.totalVisibilityId_dic[(a1,a2)] = bll

        #correct the orders of pairs in goodpair
        def correct_pairorder(pair):
            ####try:
                ####self.totalVisibilityId.tolist().index([pair[0],pair[1]])
                ####return True
            ####except:
                ####try:
                    ####self.totalVisibilityId.tolist().index([pair[1], pair[0]])
                    ####return False
                ####except:
                    ####return None
            if (pair[0], pair[1]) in self.totalVisibilityId_dic:
                return True
            elif (pair[1], pair[0]) in self.totalVisibilityId_dic:
                return False
            else:
                return None
        if verbose:
            timer.tick('c')
        #exclude pairs that are not in totalVisibilityId
        temp = []
        for p in goodpairs:
            cond = correct_pairorder([subsetant[p[0]],subsetant[p[1]]])
            if cond == True:
                temp.append(p)
            if cond == False:
                #print "correcting"
                temp.append(p[::-1])
        goodpairs = temp

        #goodpairs = [correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) for p in goodpairs if (correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) != None and correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) == True)]  #correct_pairorder([subsetant[p[0]],subsetant[p[1]]])
        nBaseline=len(goodpairs)
        if verbose:
            timer.tick('c')
        #from a pair of good antenna index to baseline index
        subsetbl = np.array([self.get_baseline([subsetant[bl[0]],subsetant[bl[1]]]) for bl in goodpairs])
        if verbose:
            timer.tick('c')
        ##################################################################################
        #bltoubl: cross bl number to ubl index
        ####def findublindex(pair,ubl=ubl):
            ####i=pair[0]
            ####j=pair[1]
            ####for k in range(len(ubl)):
                ####if dis(antloc[i]-antloc[j],ubl[k])<tolerance or dis(antloc[i]-antloc[j],-ubl[k])<tolerance:
                    ####return k
            ####print pair
            ####return "no match"
        def findublindex(pair):
            if (subsetant[pair[0]], subsetant[pair[1]]) in self.totalVisibilityUBL:
                return self.ubl2goodubl[self.totalVisibilityUBL[(subsetant[pair[0]], subsetant[pair[1]])]]
        bltoubl=[];
        for p in goodpairs:
            if p[0]!=p[1]:
                bltoubl.append(findublindex(p))
        if verbose:
            timer.tick('d')
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
                reverse.append(-1)
            elif dis(antloc[i]-antloc[j],-ubl[bltoubl[k]])<tolerance:
                reverse.append(1)
            else :
                print "something's wrong with bltoubl", crosspair[k], antloc[i]-antloc[j], bltoubl[k], ubl[bltoubl[k]]
                print i,j, subsetant[i], subsetant[j]
                print self.totalVisibilityUBL[(subsetant[i], subsetant[j])]
                exit(1)
        if verbose:
            timer.tick('e')
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
        if verbose:
            timer.tick('f')
        ###################################################
        #bl2d:  from 1d bl index to a pair of antenna numbers
        bl2d=[]
        for pair in goodpairs:
            bl2d.append(pair)#(pair[::-1])
        bl2d=np.array(bl2d)
        if verbose:
            timer.tick('g')
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
        if verbose:
            timer.tick('h')
        ####################################################################################
        #ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
        countdict={}
        for bl in bltoubl:
            countdict[bl]=[]

        for i in range(len(crosspair)):
            ant1=crosspair[i][0]
            ant2=crosspair[i][1]
            countdict[bltoubl[i]].append([ant1,ant2,i])  #([ant1,ant2,i])

        # XXX clean this up
        ublindex=[]
        for i in range(nUBL):
            ublindex.append(countdict[i])
        #turn each list in ublindex into np array
        for i in range(len(ublindex)):
            ublindex[i]=np.array(ublindex[i])
        #ublindex=np.array(ublindex)
        ublindex=np.concatenate(ublindex).astype(np.int32)
        if verbose:
            timer.tick('i')
        ###############################################################################
        #bl1dmatrix: a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
                #I suppose 99999 for bad and auto baselines?
        bl1dmatrix=(2**31-1)*np.ones([nAntenna,nAntenna],dtype='int32')
        for i in range(len(crosspair)):
            bl1dmatrix[crosspair[i][1]][crosspair[i][0]]=i
            bl1dmatrix[crosspair[i][0]][crosspair[i][1]]=i
        if verbose:
            timer.tick('j')
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

        m1=-a.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
        m2=d.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
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
            B[i][crosspair[i][0]]=reverse[i]*-1   #1
            B[i][crosspair[i][1]]=reverse[i]*1  #-1
            B[i][nAntenna+bltoubl[i]]=1
        B=sps.csr_matrix(B)
        if verbose:
            timer.tick('k')
        ###########################################################################
        #create info dictionary
        #info={}
        # XXX could interleave info assignment with variables above
        info = RedundantInfo()
        info['nAntenna']=nAntenna
        info['nUBL']=nUBL
        info['nBaseline']=nBaseline
        info['subsetant'] = np.array(subsetant,dtype=np.int32)
        info['antloc'] = np.array(antloc, dtype=np.float32)
        info['subsetbl'] = np.array(subsetbl, dtype=np.int32)
        info['ubl'] = np.array(ubl, dtype=np.float32)
        info['bltoubl'] = np.array(bltoubl, dtype=np.int32)
        info['reversed'] = np.array(reverse, dtype=np.int32)
        info['reversedauto'] = np.array(reversedauto, dtype=np.int32)
        info['autoindex'] = np.array(autoindex, dtype=np.int32)
        info['crossindex'] = np.array(crossindex, dtype=np.int32)
        #info['ncross']=ncross
        info['bl2d'] = np.array(bl2d, dtype=np.int32)
        info['ublcount'] = np.array(ublcount, dtype=np.int32)
        info['ublindex'] = ublindex
        info['bl1dmatrix'] = np.array(bl1dmatrix, dtype=np.int32)
        info['degenM'] = np.array(degenM, dtype=np.float32)
        info['At'] = A.T
        info['Bt'] = B.T
        if verbose:
            timer.tick('l')
        info.update()
        '''
        with warnings.catch_warnings():
                warnings.filterwarnings("ignore",category=DeprecationWarning)
                info['At'] = A.transpose()
                info['Bt'] = B.transpose()
                if verbose:
                    timer.tick('m')
                info['AtAi'] = la.pinv(info['At'].dot(info['A']).todense(), cond = 10**(-6))#(AtA)^-1
                info['BtBi'] = la.pinv(info['Bt'].dot(info['B']).todense(), cond = 10**(-6))#(BtB)^-1
                #if verbose:
                    #timer.tick('m')
                #info['AtAiAt'] = info['AtAi'].dot(info['At'].todense())#(AtA)^-1At
                #info['BtBiBt'] = info['BtBi'].dot(info['Bt'].todense())#(BtB)^-1Bt
                #if verbose:
                    #timer.tick('m')
                #info['PA'] = info['A'].dot(info['AtAiAt'])#A(AtA)^-1At
                #info['PB'] = info['B'].dot(info['BtBiBt'])#B(BtB)^-1Bt
                #if verbose:
                    #timer.tick('m')
                #info['ImPA'] = sps.identity(ncross) - info['PA']#I-PA
                #info['ImPB'] = sps.identity(ncross) - info['PB']#I-PB
        info['totalVisibilityId'] = self.totalVisibilityId
        '''
        if verbose:
            timer.tick('m')
        self.Info = info
        if verbose:
            timer.tick('n')


    def get_baseline(self,pair):
        '''inverse function of totalVisibilityId, calculate the baseline index from 
        the antenna pair. It allows flipping of a1 and a2, will return same result'''
        if not (type(pair) == list or type(pair) == np.ndarray or type(pair) == tuple):
            raise Exception("input needs to be a list of two numbers")
            return
        elif len(np.array(pair)) != 2:
            raise Exception("input needs to be a list of two numbers")
            return
        elif type(pair[0]) == str or type(pair[0]) == np.string_:
            raise Exception("input needs to be number not string")
            return
        ####try:
            ####return self.totalVisibilityId.tolist().index([pair[0],pair[1]])
        ####except:
            ####try:
                ####return self.totalVisibilityId.tolist().index([pair[1], pair[0]])
            ####except:
                #####raise Exception("Error: antenna pair %s not found in self.totalVisibilityId."%pair)
                ####return None
        if self.totalVisibilityId_dic is None:
            self.totalVisibilityId_dic = {}
            for bll, (a1, a2) in enumerate(self.totalVisibilityId):
                self.totalVisibilityId_dic[(a1,a2)] = bll
        if (pair[0],pair[1]) in self.totalVisibilityId_dic:
            return self.totalVisibilityId_dic[(pair[0],pair[1])]
        elif (pair[1],pair[0]) in self.totalVisibilityId_dic:
            return self.totalVisibilityId_dic[(pair[1],pair[0])]
        else:
            return None


    def compute_UBL_old2(self,tolerance = 0.1):
        '''compute_UBL returns the average of all baselines in that ubl group'''
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
            #for j in range(i+1,len(antloc)):    #(this gives the same redundant info as the correct info saved in test)
            for j in range(i+1):
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


    def compute_UBL_old(self,tolerance = 0.1):
        '''compute_UBL returns the average of all baselines in that ubl group'''
        #check if the tolerance is not a string
        if type(tolerance) == str:
            raise Exception("tolerance needs to be number not string")
            return
        ubllist = np.array([np.array([np.array([0,0,0]),1])]);
        for pair in self.totalVisibilityId:
            if pair[0] not in self.badAntenna and pair[1] not in self.badAntenna:
                [i,j] = pair
                bool = True
                for k in range(len(ubllist)):
                    if  la.norm(self.antennaLocation[i]-self.antennaLocation[j]-ubllist[k][0])<tolerance:
                        n=ubllist[k][1]
                        ubllist[k][0]=1/(n+1.0)*(n*ubllist[k][0]+self.antennaLocation[i]-self.antennaLocation[j])
                        ubllist[k][1]+=1
                        bool = False
                    elif  la.norm(self.antennaLocation[i]-self.antennaLocation[j]+ubllist[k][0])<tolerance:
                        n=ubllist[k][1]
                        ubllist[k][0]=1/(n+1.0)*(n*ubllist[k][0]-(self.antennaLocation[i]-self.antennaLocation[j]))
                        ubllist[k][1]+=1
                        bool = False
                if bool :
                    ubllist = np.append(ubllist,np.array([np.array([self.antennaLocation[j]-self.antennaLocation[i],1])]),axis=0)
        ubllist = np.delete(ubllist,0,0)
        ublall=[]
        for ubl in ubllist:
            ublall.append(ubl[0])
        ublall=np.array(ublall)
        return ublall

    def compute_UBL(self,tolerance = 0.1):
        '''XXX DOCSTRING'''
        if tolerance == 0:
            tolerance = np.min(np.linalg.norm(np.array(self.antennaLocation) - self.antennaLocation[0], axis = 1)) / 1.e6
        ubl = {}
        for bl, (a1, a2) in enumerate(self.totalVisibilityId):
            if a1 != a2 and a1 not in self.badAntenna and a2 not in self.badAntenna:
                loc_tuple = tuple(np.round((self.antennaLocation[a2] - self.antennaLocation[a1]) / float(tolerance)) * tolerance)
                neg_loc_tuple = tuple(np.round((self.antennaLocation[a1] - self.antennaLocation[a2]) / float(tolerance)) * tolerance)
                if loc_tuple in ubl:
                    ubl[loc_tuple].add(bl + 1)
                elif neg_loc_tuple in ubl:
                    ubl[neg_loc_tuple].add(- bl - 1)
                else:
                    if loc_tuple[0] >= 0:
                        ubl[loc_tuple] = set([bl + 1])
                    else:
                        ubl[neg_loc_tuple] = set([-bl - 1])

        #calculate actual average of the gridded baseline vectors to get an accurate representation of the ubl vector
        ubl_vec = np.zeros((len(ubl), 3))
        self.totalVisibilityUBL = {}

        ublcount = np.zeros(len(ubl))
        for u, grid_ubl_vec in enumerate(ubl):
            for bl in ubl[grid_ubl_vec]:
                assert bl != 0
                a1, a2 = self.totalVisibilityId[abs(bl) - 1]
                if bl > 0:
                    ubl_vec[u] = ubl_vec[u] + self.antennaLocation[a2] - self.antennaLocation[a1]
                else:
                    ubl_vec[u] = ubl_vec[u] + self.antennaLocation[a1] - self.antennaLocation[a2]
                self.totalVisibilityUBL[(a1, a2)] = u
            ublcount[u] = len(ubl[grid_ubl_vec])
            ubl_vec[u] = ubl_vec[u] / ublcount[u]

        reorder = (ubl_vec[:,1]*1e9 + ubl_vec[:,0]).argsort()
        rereorder = reorder.argsort()
        for key in self.totalVisibilityUBL:
            self.totalVisibilityUBL[key] = rereorder[self.totalVisibilityUBL[key]]
        ubl_vec = ubl_vec[reorder]

        #now I need to deal with the fact that no matter how coarse my grid is, it's possible for a single group of ubl to fall into two adjacent grids. So I'm going to check if any of the final ubl vectors are seperated by less than tolerance. If so, merge them
        ublmap = {}
        for u1 in range(len(ubl_vec)):
            for u2 in range(u1):
                if la.norm(ubl_vec[u2] - ubl_vec[u1]) < tolerance or la.norm(ubl_vec[u2] + ubl_vec[u1]) < tolerance:
                    ublmap[u1] = u2
                    ubl_vec[u2] = (ubl_vec[u1] * ublcount[u1] + ubl_vec[u2] * ublcount[u2]) / (ublcount[u1] + ublcount[u2])
                    break
            ublmap[u1] = u1

        merged_ubl_vec = []
        for u in range(len(ubl_vec)):
            if ublmap[u] == u:
                merged_ubl_vec.append(ubl_vec[u])
                ublmap[u] = len(merged_ubl_vec) - 1
            else:
                ublmap[u] = ublmap[ublmap[u]]
        merged_ubl_vec = np.array(merged_ubl_vec)

        for key in self.totalVisibilityUBL:
            self.totalVisibilityUBL[key] = ublmap[self.totalVisibilityUBL[key]]
        return ubl_vec


    def get_ublindex(self,antpair):
        '''need to do compute_redundantinfo first for this function to work (needs 'bl1dmatrix')
        input the antenna pair(as a list of two numbers), return the corresponding ubl index'''
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


    def get_reversed(self,antpair):
        '''need to do compute_redundantinfo first
        input the antenna pair, return -1 if it is a reversed baseline and 1 if it is not reversed'''
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

##########################Sub-class#############################
# XXX application to PAPER should be in another file
class RedundantCalibrator_PAPER(RedundantCalibrator):
    '''XXX DOCSTRING'''
    def __init__(self, aa):
        nTotalAnt = len(aa)
        RedundantCalibrator.__init__(self, nTotalAnt)
        self.aa = aa
        self.antennaLocationAtom = np.zeros((self.nTotalAnt,3), dtype='int32')
        for i in range(len(self.aa.ant_layout)):
            for j in range(len(self.aa.ant_layout[0])):
                self.antennaLocationAtom[self.aa.ant_layout[i][j]] = np.array([i, j, 0])

        self.preciseAntennaLocation = .299792458 * np.array([ant.pos for ant in self.aa])
        self._goodAntenna = self.aa.ant_layout.flatten()
        self._goodAntenna.sort()
        self.badAntenna = []
        self.badUBLpair = []
        for i in range(nTotalAnt):
            if i not in self._goodAntenna:
                self.badAntenna.append(i)
        self.antennaLocationAtom = self.antennaLocationAtom - np.mean(self.antennaLocationAtom[self._goodAntenna], axis = 0).astype('int32')

        ##self.antennaLocation = np.copy(self.antennaLocationAtom) #* [4, 30, 0]
        ####fit for idealized antloc
        A = np.array([list(a) + [1] for a in self.antennaLocationAtom[self._goodAntenna]])
        self.antennaLocation = np.zeros_like(self.antennaLocationAtom).astype('float64')
        self.antennaLocation[self._goodAntenna] = self.antennaLocationAtom[self._goodAntenna].dot(la.pinv(A.transpose().dot(A)).dot(A.transpose().dot(self.preciseAntennaLocation[self._goodAntenna]))[:3])##The overall constant is so large that it screws all the matrix inversion up. so im not including the over all 1e8 level shift
        self.antennaLocation[self._goodAntenna, ::2] = self.antennaLocation[self._goodAntenna, ::2].dot(np.array([[np.cos(PI/2+aa.lat), np.sin(PI/2+aa.lat)],[-np.sin(PI/2+aa.lat), np.cos(PI/2+aa.lat)]]).transpose())###rotate into local coordinates

class RedundantCalibrator_X5(RedundantCalibrator):
    def __init__(self, antennaLocation):
        nant = len(antennaLocation)
        RedundantCalibrator.__init__(self, nant)
        self.antennaLocation = antennaLocation
        self.totalVisibilityId = np.concatenate([[[i,j] for j in range(i, nant)] for i in range(nant)])

        self.badAntenna = range(16) + range(56,60) + [16,19,50]

# XXX utility function belongs in another file
def read_ndarray(path, shape, dtype, ranges):
    '''read middle part of binary file of shape and dtype specified by ranges of the first dimension. ranges is [inclusive, exclusive)'''
    if not os.path.isfile(path):
        raise IOError(path + 'doesnt exist.')
    if len(ranges) != 2 or ranges[0] < 0 or ranges[0] >= ranges[1] or ranges[1] > shape[0]:
        raise ValueError("%s is not a vlid range."%ranges)
    nbytes = np.dtype(dtype).itemsize
    higher_dim_chunks = 1 # product of higher dimensions. if array is (2,3,4,5), this is 3*4*5
    for m in shape[1:]:
        higher_dim_chunks = higher_dim_chunks * m

    #print np.fromfile(path, dtype = dtype).shape
    with open(path, 'r') as f:
        f.seek(higher_dim_chunks * nbytes * ranges[0])
        #print higher_dim_chunks * nbytes * ranges[0]
        result = np.fromfile(f, dtype = dtype, count = (ranges[1] - ranges[0]) * higher_dim_chunks)
    new_shape = np.array(shape)
    new_shape[0] = (ranges[1] - ranges[0])
    #print result.shape,tuple(new_shape)
    return result.reshape(tuple(new_shape))

# XXX utility function belongs in another file
def write_ndarray(path, shape, dtype, ranges, data, check = True, max_retry = 3, task = 'unkown'):
    '''write middle part of binary file of shape and dtype specified by ranges of the first dimension. ranges is [inclusive, exclusive)'''
    if not os.path.isfile(path):
        raise IOError(path + 'doesnt exist.')
    if len(ranges) != 2 or ranges[0] < 0 or ranges[0] >= ranges[1] or ranges[1] > shape[0]:
        raise ValueError("%s is not a vlid range."%ranges)
    if data.dtype != dtype or data.shape[1:] != shape[1:] or data.shape[0] != ranges[1] - ranges[0]:
        raise ValueError("data shape %s cannot be fit into data file shape %s."%(data.shape, shape))
    nbytes = np.dtype(dtype).itemsize
    higher_dim_chunks = 1 # product of higher dimensions. if array is (2,3,4,5), this is 3*4*5
    for m in shape[1:]:
        higher_dim_chunks = higher_dim_chunks * m
    with open(path, 'r+') as f:
        f.seek(higher_dim_chunks * nbytes * ranges[0])
        data.tofile(f)
    if check:
        tries = 0
        while not (data == read_ndarray(path, shape, dtype, ranges)).all() and tries < max_retry:

            time.sleep(1)
            tries = tries + 1
            with open(path, 'r+') as f:
                f.seek(higher_dim_chunks * nbytes * ranges[0])
                data.tofile(f)
        if not (data == read_ndarray(path, shape, dtype, ranges)).all():
            raise IOError("write_ndarray failed on %s with shape %s between %s with task %s."%(path, shape, ranges, task))
    return

def load_omnichisq(path):
    '''XXX DOCSTRING'''
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise IOError("Path %s does not exist."%path)

    omnichisq = np.fromfile(path, dtype = 'float32')
    NF = int(omnichisq[2])
    omnichisq.shape = (len(omnichisq) / (NF + 3), (NF + 3))
    return omnichisq

def load_omnigain(path, info=None):
    '''XXX DOCSTRING'''
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise IOError("Path %s does not exist."%path)
    if info is None:
        info = path.replace('.omnigain', '.binfo')
    if type(info) == type('a'):
        info = read_redundantinfo(info)


    omnigain = np.fromfile(path, dtype = 'float32')
    omnigain.shape = (omnigain.shape[0] / (info['nAntenna']) / (2 + 1 + 1 + 2 * int(omnigain[3])), info['nAntenna'], 2 + 1 + 1 + 2 * int(omnigain[3]))
    return omnigain

def load_omnifit(path, info=None):
    '''XXX DOCSTRING'''
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise IOError("Path %s does not exist."%path)
    if info is None:
        info = path.replace('.omnifit', '.binfo')
    if type(info) == type('a'):
        info = read_redundantinfo(info)


    omnifit = np.fromfile(path, dtype = 'float32')
    omnifit.shape = (omnifit.shape[0] / (info['nUBL']) / (2 + 3 + 1 + 2 * int(omnifit[5])), info['nUBL'], 2 + 3 + 1 + 2 * int(omnifit[3]))
    return omnifit

def get_omnitime(omnistuff):
    '''XXX DOCSTRING'''
    if len(omnistuff.shape) == 2:
        return np.array([struct.unpack('d', struct.pack('ff', *(pair.tolist())))[0] for pair in omnistuff[:, :2]])
    elif len(omnistuff.shape) == 3:
        return np.array([struct.unpack('d', struct.pack('ff', *(pair.tolist())))[0] for pair in omnistuff[:, 0, :2]])
    else:
        raise ValueError('get_omnitime does not know how to deal with array of shape %s.'%omnistuff.shape)

# XXX utility function belongs in another file
def omniview(data_in, info, plotrange = None, oppath = None, suppress = False, title = '', plot_single_ubl = False, plot_3 = False, plot_1 = -1):
    '''plot_3: only plot the 3 most redundant ones. plot_1: counting start from 0 the most redundant baseline'''
    import matplotlib.pyplot as plt
    data = np.array(data_in)
    try:#in case info is Info class
        info = info.get_info()
    except:
        pass
    if plot_3 and info['nUBL'] < 3:
        plot_3 = False

    colors=[]
    colorgrid = int(math.ceil((info['nUBL']/12.+1)**.34))
    for red in range(colorgrid):
        for green in range(colorgrid):
            for blue in range(colorgrid):
                #print red, green, blue
                colors += [(np.array([red, green, blue]).astype('float')/(colorgrid - 1)).tolist()]
    #colors.remove([0,0,0])
    colors.remove([1,1,1])

    if plot_3:
        select_ubl_index = np.argsort(info['ublcount'])[-3:]
    elif plot_1 >= 0:
        select_ubl_index = np.argsort(info['ublcount'])[::-1][plot_1:plot_1+1]


    if len(data.shape) == 1 or data.shape[0] == 1:
        ds = [data[info['subsetbl']][info['crossindex']]]
        fig, axes = plt.subplots(nrows=1, ncols=1, sharey=True, sharex=True)
        axes = [axes]
    else:
        ds = data[:, info['subsetbl'][info['crossindex']]]
        fig, axes = plt.subplots(nrows=1, ncols=len(ds), sharey=True, sharex=True)

    outputdata = []
    for i in range(len(ds)):
        outputdata = outputdata + [[]]
        d = ds[i]
        ax = axes[i]
        if plotrange is None:
            plotrange = 1.2*np.nanmax(np.abs(d))

        ubl = 0
        for marker in ["o", "v", "^", "<", ">", "8", "s", "p", "h", (6,1,0), (8,1,0), "d"]:
            for color in colors:
                #print info['ublindex'][ubl][:,2]
                #print marker, color
                if (plot_single_ubl or len(info['ublindex'][ubl]) > 1) and (not (plot_3 or plot_1 >= 0) or ubl in select_ubl_index):
                    if plot_3 or plot_1 >= 0:
                        color = [[1,0,0],[0,1,0],[0,0,1]][select_ubl_index.tolist().index(ubl)]
                    ax.scatter(np.real(d[np.array(info['ublindex'][ubl][:,2]).astype('int')]),np.imag(d[np.array(info['ublindex'][ubl][:,2]).astype('int')])*info['reversed'][np.array(info['ublindex'][ubl][:,2]).astype('int')], marker=marker, color=color)
                    outputdata[i] = outputdata[i] + [(np.real(d[np.array(info['ublindex'][ubl][:,2]).astype('int')]) + 1.j * np.imag(d[np.array(info['ublindex'][ubl][:,2]).astype('int')])*info['reversed'][np.array(info['ublindex'][ubl][:,2]).astype('int')], marker, color, info['ubl'][ubl])]

                ubl += 1
                if ubl == info['nUBL']:
                    #if i == 1:
                        #ax.text(-(len(ds)-1 + 0.7)*plotrange, -0.7*plotrange, "#Ant:%i\n#UBL:%i"%(info['nAntenna'],info['nUBL']),bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))
                    ax.set_title(title + "\nGood Antenna count: %i\nUBL count: %i"%(info['nAntenna'],info['nUBL']))
                    ax.grid(True)
                    ax.set(adjustable='datalim', aspect=1)
                    ax.set_xlabel('Real')
                    ax.set_ylabel('Imag')
                    break
            if ubl == info['nUBL']:
                break
    plt.axis([-plotrange, plotrange, -plotrange, plotrange])
    if oppath is not None:
        plt.savefig(oppath, bbox_inches='tight')
    if not suppress:
        plt.show()
    else:
        plt.close()
    return outputdata

# XXX utility function belongs in another file
def lin_depend(v1, v2, tol = 0):
    '''whether v1 and v2 are linearly dependent'''
    if len(v1) != len(v2):
        raise Exception("Length mismatch %i vs %i."%(len(v1), len(v2)))
    if la.norm(v1) == 0:
        return True
    return la.norm(np.dot(v1, v2)/np.dot(v1, v1) * v1 - v2) <= tol

# XXX utility function belongs in another file
def _f(rawcal_ubl=[], verbose=False):
    '''run this function twice in a row and its christmas'''
    if verbose and rawcal_ubl != []:
        print "Starting ubl:", rawcal_ubl
    if rawcal_ubl == []:
        rawcal_ubl += [2,3]
    if verbose:
        print "ubl:", rawcal_ubl

def find_solution_path(info, input_rawcal_ubl=[], tol = 0.0, verbose=False):
    '''return (intialantenna, solution_path) for raw calibration. solution path
    contains a list of [(a0, a1, crossubl), a] = [(ublindex entry), (which ant is
    solved, 0 or 1)]. When raw calibrating, initialize calpar to have [0] at
    initial antenna, then simply iterate along the solution_path, use crossubl and
    a0 or a1 specified by a to solve for the other a1 or a0 and append it to
    calpar. Afterwards, use mean angle on calpars'''
    ###select 2 ubl for calibration
    rawcal_ubl = list(input_rawcal_ubl)
    if verbose and rawcal_ubl != []:
        print "Starting ubl:", rawcal_ubl
    if rawcal_ubl == []:
        ublcnt_tmp = info['ublcount'].astype('float')
        rawcal_ubl += [np.argmax(ublcnt_tmp)]
        if verbose:
            print "Picking %s with redundancy %i as first ubl"%(info['ubl'][rawcal_ubl[-1]], ublcnt_tmp[rawcal_ubl[-1]])
        ublcnt_tmp[rawcal_ubl[-1]] = np.nan
        rawcal_ubl += [np.nanargmax(ublcnt_tmp)]
        if verbose:
            print "Picking %s with redundancy %i as second ubl"%(info['ubl'][rawcal_ubl[-1]], ublcnt_tmp[rawcal_ubl[-1]])
        ublcnt_tmp[rawcal_ubl[-1]] = np.nan
        #while np.allclose(info['ubl'][rawcal_ubl[0]]/(la.norm(info['ubl'][rawcal_ubl[0]])/la.norm(info['ubl'][rawcal_ubl[1]])), info['ubl'][rawcal_ubl[1]]) or np.allclose(info['ubl'][rawcal_ubl[0]]/(la.norm(info['ubl'][rawcal_ubl[0]])/la.norm(info['ubl'][rawcal_ubl[1]])), -info['ubl'][rawcal_ubl[1]]):
        while lin_depend(info['ubl'][rawcal_ubl[0]], info['ubl'][rawcal_ubl[1]], tol=tol):
            if verbose:
                print info['ubl'][rawcal_ubl[0]], "and", info['ubl'][rawcal_ubl[1]], "are linearly dependent."
            try:
                rawcal_ubl[1] = np.nanargmax(ublcnt_tmp)
                if verbose:
                    print "Picking %s with redundancy %i as second ubl"%(info['ubl'][rawcal_ubl[-1]], ublcnt_tmp[rawcal_ubl[-1]])
            except:
                raise Exception("Cannot find two unique baselines that are linearly independent!")
            ublcnt_tmp[rawcal_ubl[-1]] = np.nan
    if verbose:
        print "ubl:", info['ubl'][rawcal_ubl[0]], info['ubl'][rawcal_ubl[1]]

    if info['ublcount'][rawcal_ubl[0]] + info['ublcount'][rawcal_ubl[1]] <= info['nAntenna'] + 2:
        raise Exception('Array not redundant enough! Two most redundant baselines %s and %s have %i and %i baselines, which is not larger than 2 + %i'%(info['ubl'][rawcal_ubl[0]],info['ubl'][rawcal_ubl[1]], info['ublcount'][rawcal_ubl[0]],info['ublcount'][rawcal_ubl[1]], info['nAntenna']))

    ublindex = np.concatenate((np.array(info['ublindex'][rawcal_ubl[0]]).astype('int'), np.array(info['ublindex'][rawcal_ubl[1]]).astype('int')))#merge ublindex since we set both ubl phase to 0

    ###The overarching goal is to find a solution path (a sequence of unique baselines to solve) that can get multiple solutions to multiple antennas using just two sets of ubls
    solution_path = []

    antcnt = np.bincount(np.array(ublindex)[:,:2].flatten(), minlength = info['nAntenna'])#how many times each antenna appear in the two sets of ubls. at most 4. not useful if < 2
    unsolved_ant = []
    for a in range(len(antcnt)):
        if antcnt[a] == 0:
            unsolved_ant.append(a)
    if verbose:
        print "antcnt", antcnt, "Antennas", np.array(info['subsetant'])[unsolved_ant], "not directly solvable."


    ###Status string for ubl: NoUse: none of the two ants have been solved; Solvable: at least one of the ants have solutions; Done: used to generate one antennacalpar
    ublstatus = ["NoUse" for i in ublindex]

    ###antenna calpars, a list for each antenna
    calpar = np.array([[]] * info['nAntenna']).tolist()
    ###select initial antenna
    initialant = int(np.argmax(antcnt))
    if verbose:
        print "initialant", np.array(info['subsetant'])[initialant]
    calpar[initialant].append(0)
    for i in range(len(ublstatus)):
        if initialant in ublindex[i, 0:2]:
            ublstatus[i] = "Solvable"

    ###start looping
    solvecnt = 10#number of solved baselines in each loop, 10 is an arbitrary starting point
    if verbose:
        print "new ant solved",
    while solvecnt > 0:
        solvecnt = 0
        for i in range(len(ublstatus)):
            if ublstatus[i] == "Solvable":
                solvecnt += 1
                if calpar[ublindex[i, 0]] != []:#if the first antenna is solved
                    #print ublindex[i], ublindex[i, 1], len(calpar[ublindex[i, 1]]), calpar[ublindex[i, 1]],
                    calpar[ublindex[i, 1]].append(0)#just append a dummy
                    ublstatus[i] = "Done"
                    solution_path.append([ublindex[i], 0])
                    #print len(calpar[ublindex[i, 1]]), calpar[ublindex[i, 1]]
                    if len(calpar[ublindex[i, 1]]) == 1:
                        if verbose:
                            print np.array(info['subsetant'])[ublindex[i, 1]],
                        for j in range(len(ublstatus)):
                            if (ublindex[i, 1] in ublindex[j, 0:2]) and ublstatus[j] == "NoUse":
                                ublstatus[j] = "Solvable"
                else:
                    #print ublindex[i], ublindex[i, 0], len(calpar[ublindex[i, 0]]), calpar[ublindex[i, 0]],
                    calpar[ublindex[i, 0]].append(0)#just append a dummy
                    ublstatus[i] = "Done"
                    solution_path.append([ublindex[i], 1])
                    #print len(calpar[ublindex[i, 0]]), calpar[ublindex[i, 0]]
                    if len(calpar[ublindex[i, 0]]) == 1:
                        if verbose:
                            print np.array(info['subsetant'])[ublindex[i, 0]],
                        for j in range(len(ublstatus)):
                            if (ublindex[i, 0] in ublindex[j, 0:2]) and ublstatus[j] == "NoUse":
                                ublstatus[j] = "Solvable"
    if verbose:
        print ""
        if len(solution_path) != len(ublindex):
            print "Solution path has %i entries where as total candidates in ublindex have %i. The following baselines form their isolated isaland:"%(len(solution_path), len(ublindex))
            unsolved_ubl = []
            for i in range(len(ublstatus)):
                if ublstatus[i] != "Done":
                    print np.array(info['subsetant'])[ublindex[i][0]], np.array(info['subsetant'])[ublindex[i][1]]
                    unsolved_ubl.append(ublindex[i])
            unsolved_ubl = np.array(unsolved_ubl)[:,:2].flatten()
            for a in range(info['nAntenna']):
                if a in unsolved_ubl:
                    unsolved_ant.append(a)
    ant_solved = 10
    additional_solution_path = []
    while len(unsolved_ant) > 0 and ant_solved > 0:#find a ubl that can solve these individual antennas not involved in the chosen 2 ubls. Use while because the first ant in the unsolved_ant may not be solvable on the first pass

        ant_solved = 0
        for a in unsolved_ant:
            if verbose:
                print "trying to solve for ", np.array(info['subsetant'])[a],
            ublcnt_tmp = info['ublcount'].astype('float')
            third_ubl_good = False
            tried_all_ubl = False
            while (not third_ubl_good) and (not tried_all_ubl):
                try:
                    third_ubl = np.nanargmax(ublcnt_tmp)
                    ublcnt_tmp[third_ubl] = np.nan
                except:
                    tried_all_ubl = True
                    break
                if verbose:
                    print "trying ubl ", third_ubl,
                third_ubl_good = False #assume false and start checking if this ubl 1) has this antenna 2) has another baseline whose two ants are both solved
                if (len(info['ublindex'][third_ubl]) < 2) or (a not in info['ublindex'][third_ubl]):
                    continue
                third_ubl_good1 = False
                third_ubl_good2 = False
                for a1, a2, bl in info['ublindex'][third_ubl]:
                    if (a1 not in unsolved_ant) and (a2 not in unsolved_ant):
                        third_ubl_good1 = True
                        if third_ubl_good2:
                            break
                    if ((a == a1) and (a2 not in unsolved_ant)) or ((a == a2) and (a1 not in unsolved_ant)):
                        third_ubl_good2 = True
                        if third_ubl_good1:
                            break
                third_ubl_good = (third_ubl_good1 and third_ubl_good2)
            if third_ubl_good:#figure out how to use this third ubl to solve this a
                if verbose:
                    print "picked ubl", info['ubl'][third_ubl], "to solve for ant", np.array(info['subsetant'])[a]
                get_ubl_fit = []#a recipe for how to get the ublfit and solvefor the unsolved antenna
                for a1, a2, bl in info['ublindex'][third_ubl].astype('int'):
                    if (a1 not in unsolved_ant) and (a2 not in unsolved_ant):
                        get_ubl_fit.append([a1, a2, bl, info['reversed'][bl]])
                for a1, a2, bl in info['ublindex'][third_ubl].astype('int'):
                    if (a1 not in unsolved_ant) and (a2 == a):
                        get_ubl_fit.append([a1, a2, bl, info['reversed'][bl], 0])
                        break
                    if (a2 not in unsolved_ant) and (a1 == a):
                        get_ubl_fit.append([a1, a2, bl, info['reversed'][bl], 1])
                        break
                additional_solution_path.append(get_ubl_fit)
                ant_solved += 1
                unsolved_ant.remove(a)

    #remove the effect of enforcing the two baselines to be 0, rather, set the first two linearly independent antennas w.r.t initant to be 0
    # find two antennas:
    a1 = initialant
    a2 = initialant
    for index in info['ublindex'][rawcal_ubl[0]]:
        if index[0] == initialant:
            a1 = int(index[1])
            break
        elif index[1] == initialant:
            a1 = int(index[0])
            break
    bl1 = np.array(info['antloc'][a1]) - info['antloc'][initialant]
    for index in info['ublindex'][rawcal_ubl[1]]:
        if index[0] == initialant:
            a2 = int(index[1])
            break
        elif index[1] == initialant:
            a2 = int(index[0])
            break
    bl2 = np.array(info['antloc'][a2]) - info['antloc'][initialant]
    #A = np.array([bl1[:2], bl2[:2]])
    #remove_Matrix = (np.array(info['antloc'])- info['antloc'][initialant])[:,:2].dot(la.pinv(A.transpose().dot(A)).dot(A.transpose()))
    A = np.array([bl1, bl2])
    remove_Matrix = (np.array(info['antloc'])- info['antloc'][initialant]).dot(la.pinv(A.transpose().dot(A)).dot(A.transpose()))
    degeneracy_remove = [a1, a2, remove_Matrix]
    if verbose:
        print "Degeneracy: a1 = %i, a2 = %i"%(info['subsetant'][a1], info['subsetant'][a2])
    return initialant, solution_path, additional_solution_path, degeneracy_remove, (unsolved_ant == [])

def meanAngle(a, weights = None, axis = -1):
    return np.angle(np.average(np.exp(1.j*np.array(a)), weights = weights, axis = axis))

def medianAngle(a, axis = -1):
    return np.angle(nanmedian(np.cos(a), axis = axis) + 1.j * nanmedian(np.sin(a), axis = axis))

#def _medianAngle(data, result, axis = -1):
    #result_shape = collapse_shape(data.shape, axis)

    #np_result = np.frombuffer(result, dtype='float32')
    #np_result.shape = tuple(result_shape)
    #np_result[:] = medianAngle(data, axis = axis).reshape(result_shape)
    #return

# XXX utility function belongs in another file
def collapse_shape(shape, axis):
    '''XXX DOCSTRING'''
    if axis == 0 or axis == -len(shape):
        return tuple(list(shape)[1:])
    elif axis == -1 or axis == len(shape) - 1:
        return tuple(list(shape)[:-1])
    else:
        return tuple(list(shape)[:axis] + list(shape)[axis+1:])

###curerntly suffering from slow initialization which is probably due to copying data into shared array. worth further investigation.
#def medianAngle_multithread(data, axis = -1, nthread = None, verbose = False):
    #if axis < 0:
        #axis = data.ndim + axis
    #parallel_axis2 = np.argmax(collapse_shape(data.shape, axis))#the axis after averaging
    #if parallel_axis2 >= axis:
        #parallel_axis1 = parallel_axis2 + 1
    #else:
        #parallel_axis1 = parallel_axis2
    #parallel_axis_len = data.shape[parallel_axis1]
    #if nthread is None:
        #nthread = min(mp.cpu_count() - 1, parallel_axis_len)
    #nthread = min(nthread, parallel_axis_len)
    #if nthread < 2 or data.ndim == 1:
        #return medianAngle(data, axis=axis)
    #else:
        #results = {}
        #np_results = {}

        #threads = {}
        #fchunk = {}
        #chunk = parallel_axis_len / int(nthread)
        #excess = parallel_axis_len % int(nthread)
        #kwarg = {"axis": axis}

####set up threads
        #for i in range(nthread):
            #if excess == 0:
                #fchunk[i] = (i * chunk, min((1 + i) * chunk, parallel_axis_len),)
            #elif i < excess:
                #fchunk[i] = (i * (chunk+1), min((1 + i) * (chunk+1), parallel_axis_len),)
            #else:
                #fchunk[i] = (fchunk[i-1][1], min(fchunk[i-1][1] + chunk, parallel_axis_len),)

            #result_shape = list(collapse_shape(data.shape, axis))
            #result_shape[parallel_axis2] = fchunk[i][1] - fchunk[i][0]

            #results[i] = mp.RawArray('f', np.prod(result_shape))
            #np_results[i] = np.frombuffer(results[i], dtype='float32')
            #np_results[i].shape = tuple(result_shape)
            #def _slice(a):
                #return a[fchunk[i][0]:fchunk[i][1]]
            #threads[i] = mp.Process(target = _medianAngle, args = (np.apply_along_axis(_slice, parallel_axis1, data), results[i]), kwargs=kwarg)

####start processing
        #if verbose:
            #print "Starting medianAngle Process",
            #sys.stdout.flush()
        #for i in range(nthread):
            #if verbose:
                #print "#%i"%i,
                #sys.stdout.flush()
            #threads[i].start()
        #if verbose:
            #print "Finished Process",
        #for i in range(nthread):
            #threads[i].join()
            #if verbose:
                #print "#%i"%i,
        #if verbose:
            #print ""
            #sys.stdout.flush()
        #return np.concatenate([np_results[i] for i in range(nthread)],axis=parallel_axis2)

def raw_calibrate(data, info, initant, solution_path, additional_solution_path, degeneracy_remove):
    '''XXX DOCSTRING'''
    result = np.ones(int(math.floor((len(data)*2.)**.5)), dtype='complex64')
    calpar = np.array([[]]*info['nAntenna']).tolist()
    calpar[initant] = [0]
    d=np.angle(data[info['subsetbl']][info['crossindex']])
    #d=np.angle(omni.apply_calpar(data, result, visibilityID)[info['subsetbl']][info['crossindex']])
    for ublindex, a in solution_path:
        calpar[ublindex[1-a]].append(calpar[ublindex[a]][0] + ((a-0.5)/-.5)*d[ublindex[2]])
    for i in range(len(calpar)):
        if len(calpar[i]) > 0:
            calpar[i] = [meanAngle(calpar[i])]
        else:
            calpar[i] = [0]

    #now deal with additional_solution_path which deal with antennas that are not included in the 2 ubls picked to be 0
    for solution in additional_solution_path:
        ubl_phases = np.array([s[-1]*(calpar[s[0]][0]-calpar[s[1]][0]+d[s[2]]) for s in solution[:-1]])
        ubl_phase = medianAngle(ubl_phases)
        #print np.angle(np.exp(1.j*ubl_phases))
        #print calpar[solution[-1][1-solution[-1][-1]]]
        calpar[solution[-1][1-solution[-1][-1]]] = [calpar[solution[-1][solution[-1][-1]]][0] + ((solution[-1][-1]-0.5)/-.5)*(d[solution[-1][2]] - ubl_phase * (solution[-1][-2]))]
        #print solution[-1]
        #print calpar[solution[-1][solution[-1][-1]]][0], ((solution[-1][-1]-0.5)/-.5),d[solution[-1][2]] , ubl_phase * (solution[-1][-2]), calpar[solution[-1][1-solution[-1][-1]]]

    calpar = (np.array(calpar).flatten() + np.pi) % (2 * np.pi) - np.pi
    #remove the effect of enforcing the two baselines to be 0, rather, set the first two linearly independent antennas w.r.t initant to be 0

    calpar = calpar - degeneracy_remove[2].dot([calpar[degeneracy_remove[0]],calpar[degeneracy_remove[1]]])

    result[info['subsetant']] = np.exp(1.j*calpar)# * result[info['subsetant']]
    return result

# XXX utility class belongs in another file
class InverseCholeskyMatrix:
    '''for a positive definite matrix, Cholesky decomposition is M = L.Lt, where L
    lower triangular. This decomposition helps computing inv(M).v faster, by
    avoiding calculating inv(M). Once we have L, the product is simply
    inv(Lt).inv(L).v, and inverse of triangular matrices multiplying a vector is
    fast. sla.solve_triangular(M, v) = inv(M).v'''
    def __init__(self, matrix):
        if type(matrix).__module__ != np.__name__ or len(matrix.shape) != 2:
            raise TypeError("matrix must be a 2D numpy array");
        try:
            self.L = la.cholesky(matrix)#L.dot(L.conjugate().transpose()) = matrix, L lower triangular
            self.Lt = self.L.conjugate().transpose()
            #print la.norm(self.L.dot(self.Lt)-matrix)/la.norm(matrix)
        except:
            raise TypeError("cholesky failed. matrix is not positive definite.")

    @classmethod
    def fromfile(cls, filename, n, dtype):
        if not os.path.isfile(filename):
            raise IOError("%s file not found!"%filename)
        matrix = cls(np.array([[1,0],[0,1]]))
        try:
            matrix.L = np.fromfile(filename, dtype=dtype).reshape((n,n))#L.dot(L.conjugate().transpose()) = matrix, L lower triangular
            matrix.Lt = matrix.L.conjugate().transpose()
            #print la.norm(self.L.dot(self.Lt)-matrix)/la.norm(matrix)
        except:
            raise TypeError("cholesky import failed. matrix is not %i by %i with dtype=%s."%(n, n, dtype))
        return matrix

    def dotv(self, vector):
        try:
            return la.solve_triangular(self.Lt, la.solve_triangular(self.L, vector, lower=True), lower=False)
        except:
            return np.empty_like(vector)+np.nan

    def dotM(self, matrix):
        return np.array([self.dotv(v) for v in matrix.transpose()]).transpose()

    def astype(self, t):
        self.L = self.L.astype(t)
        self.Lt = self.Lt.astype(t)
        return self

    def tofile(self, filename, overwrite = False):
        if os.path.isfile(filename) and not overwrite:
            raise IOError("%s file exists!"%filename)
        self.L.tofile(filename)

# XXX utility function belongs elsewhere
def solve_slope(A_in, b_in, tol, niter=30, step=1, verbose=False):
    '''solve for the solution vector x such that mod(A.x, 2pi) = b, 
    where the values range from -p to p. solution will be seeked 
    on the first axis of b'''
    timer = Timer()
    p = np.pi
    A = np.array(A_in)
    b = np.array(b_in + p) % (2*p) - p
    if A.ndim != 2:
        raise TypeError("A matrix must be 2 dimensional. Input A is %i dimensional."%A.ndim)
    if A.shape[0] != b.shape[0]:
        raise TypeError("A and b has shape mismatch: %s and %s."%(A.shape, b.shape))
    if A.shape[1] != 2:
        raise TypeError("A matrix's second dimension must have size of 2. %i inputted."%A.shape[1])
    if verbose:
        timer.tick("a")
    #find the shortest 2 non-parallel baselines, candidate_vecs have all combinations of vectors in a summation or subtraction. Each entry is i,j, v0,v1 where Ai+Aj=(v0,v1), negative j means subtraction. Identical i,j means vector itself without add or subtract
    candidate_vecs = np.zeros((len(A)**2, 4), dtype = 'float32')
    n = 0
    for i in range(len(A)):
        for j in range(len(A)):
            if i < j:
                candidate_vecs[n] = [i, j, A[i,0]+A[j,0], A[i,1]+A[j,1]]
            elif i == j:
                candidate_vecs[n] = [i, j, A[i,0], A[i,1]]
            elif i > j:
                candidate_vecs[n] = [i, -j, A[i,0]-A[j,0], A[i,1]-A[j,1]]

            n = n + 1
    if verbose:
        timer.tick("b")
    candidate_vecs = candidate_vecs[np.linalg.norm(candidate_vecs, axis=1)>tol]

    #construct coarse A that consists of the 2 shortest vecs
    coarseA = np.zeros((2,2), dtype = 'float32')
    if b.ndim > 1:
        coarseb = np.zeros(np.concatenate(([2], b.shape[1:])), dtype='float32')
    else:
        coarseb = np.zeros(2, dtype='float32')

    for n in np.argsort(np.linalg.norm(candidate_vecs[:, 2:4], axis=1)):
        v = candidate_vecs[n, 2:4]
        if la.norm(coarseA[0]) == 0:
            coarseA[0] = v
        else:
            perp_component = v - v.dot(coarseA[0])/(coarseA[0].dot(coarseA[0])) * coarseA[0]
            if la.norm(perp_component) > tol:
                coarseA[1] = v
                break
    if la.norm(coarseA[1]) == 0:
        raise Exception("Poorly constructed A matrix: cannot find a pair of orthogonal vectors")
    if verbose:
        timer.tick("c")
    #construct coarse b that contains medianAngle off all bs correponding to the 2 shortest bls
    coarseb0_candidate_indices = np.arange(len(candidate_vecs))[(np.linalg.norm(candidate_vecs[:, 2:4] - coarseA[0], axis=-1) < tol)|(np.linalg.norm(candidate_vecs[:, 2:4] + coarseA[0], axis=-1) < tol)]#stores the indices in candidate_vecs that is revelant to coarseb0
    coarseb1_candidate_indices = np.arange(len(candidate_vecs))[(np.linalg.norm(candidate_vecs[:, 2:4] - coarseA[1], axis=-1) < tol)|(np.linalg.norm(candidate_vecs[:, 2:4] + coarseA[1], axis=-1) < tol)]#stores the indices in candidate_vecs that is revelant to coarseb1
    coarseb0_candidate_shape = np.array(coarseb.shape)
    coarseb0_candidate_shape[0] = len(coarseb0_candidate_indices)
    coarseb1_candidate_shape = np.array(coarseb.shape)
    coarseb1_candidate_shape[0] = len(coarseb1_candidate_indices)
    coarseb0_candidate = np.zeros(coarseb0_candidate_shape, dtype='float32')
    coarseb1_candidate = np.zeros(coarseb1_candidate_shape, dtype='float32')

    for nn, (coarseb_candidate_indices, coarseb_candidate) in enumerate(zip([coarseb0_candidate_indices, coarseb1_candidate_indices], [coarseb0_candidate, coarseb1_candidate])):
        for n, ind in enumerate(coarseb_candidate_indices):
            i = int(candidate_vecs[ind, 0])
            j = int(candidate_vecs[ind, 1])
            v = candidate_vecs[ind, 2:4]
            if la.norm(coarseA[nn] - v) < tol:
                bsign = 1
            else:
                bsign = -1
            if i < j:
                coarseb_candidate[n] = b[i]+b[j]#(b[i]+b[j]+p)%(2*p)-p
            elif i == j:
                coarseb_candidate[n] = b[i]
            elif i > j:
                coarseb_candidate[n] = b[i]-b[abs(j)]#(b[i]-b[abs(j)]+p)%(2*p)-p
            coarseb_candidate[n] = coarseb_candidate[n] * bsign
    if verbose:
        timer.tick("d")

    coarseb[0] = medianAngle(coarseb0_candidate, axis=0)
    coarseb[1] = medianAngle(coarseb1_candidate, axis=0)
    if verbose:
        print coarseb0_candidate.shape
        timer.tick("d")
    # find coarse solutions
    try:
        icA = la.inv(coarseA)
    except:
        raise Exception("Poorly constructed coarseA matrix: %s."%(coarseA))
    try:
        #iA = InverseCholeskyMatrix(A.transpose().dot(A))
        iA = la.inv(A.transpose().dot(A))
    except:
        raise Exception("Poorly constructed A matrix: %s."%(A.transpose().dot(A)))
    if verbose:
        print iA.shape
        timer.tick("e")
    if b.ndim > 2:
        extra_shape = b.shape[1:]
        flat_extra_dim = 1
        for i in range(1, b.ndim):
            flat_extra_dim = flat_extra_dim * b.shape[i]
        coarseb.shape = (2, flat_extra_dim)
        b.shape = (len(b), flat_extra_dim)
    else:
        extra_shape = None
    if verbose:
        timer.tick("f")
    result = icA.dot(coarseb)
    if verbose:
        print coarseA
        print result
    for i in range(niter):
        result = result + step * iA.dot(A.transpose().dot((b - A.dot(result) + p)%(2*p)-p))
        if verbose:
            print result
    if extra_shape is not None:
        result.shape = tuple(np.concatenate(([2], extra_shape)))
    if verbose:
        timer.tick("g")
    return result

#def solve_slope_old(A_in, b_in, tol, niter=3, p = np.pi):#solve for the solution vector x such that mod(A.x, 2pi) = b, where the values range from -p to p. solution will be seeked on the first axis of b
    #A = np.array(A_in)
    #b = np.array(b_in + p) % (2*p) - p
    #if A.ndim != 2:
        #raise TypeError("A matrix must be 2 dimensional. Input A is %i dimensional."%A.ndim)
    #if A.shape[0] != b.shape[0]:
        #raise TypeError("A and b has shape mismatch: %s and %s."%(A.shape, b.shape))
    #if A.shape[1] != 2:
        #raise TypeError("A matrix's second dimension must have size of 2. %i inputted."%A.shape[1])

    ##find the shortest 2 non-parallel baselines, candidate_vecs have all combinations of vectors in a summation or subtraction. Each entry is i,j, v0,v1 where Ai+Aj=(v0,v1), negative j means subtraction. Identical i,j means vector itself without add or subtract
    #candidate_vecs = np.zeros((len(A)**2, 4), dtype = 'float32')
    #n = 0
    #for i in range(len(A)):
        #for j in range(len(A)):
            #if i < j:
                #candidate_vecs[n] = [i, j, A[i,0]+A[j,0], A[i,1]+A[j,1]]
            #elif i == j:
                #candidate_vecs[n] = [i, j, A[i,0], A[i,1]]
            #elif i > j:
                #candidate_vecs[n] = [i, -j, A[i,0]-A[j,0], A[i,1]-A[j,1]]

            #n = n + 1

    #candidate_vecs = candidate_vecs[np.linalg.norm(candidate_vecs, axis=1)>tol]

    ##construct coarse A that consists of the 2 shortest vecs
    #coarseA = np.zeros((2,2), dtype = 'float32')
    #if b.ndim > 1:
        #coarseb = np.zeros(np.concatenate(([2], b.shape[1:])), dtype='float32')
    #else:
        #coarseb = np.zeros(2, dtype='float32')
    #for n in np.argsort(np.linalg.norm(candidate_vecs[:, 2:4], axis=1)):
        #i = candidate_vecs[n, 0]
        #j = candidate_vecs[n, 1]
        #v = candidate_vecs[n, 2:4]
        #if la.norm(coarseA[0]) == 0:
            #coarseA[0] = v
            #if i < j:
                #coarseb[0] = (b[i]+b[j]+p)%(2*p)-p
            #elif i == j:
                #coarseb[0] = b[i]
            #elif i > j:
                #coarseb[0] = (b[i]-b[abs(j)]+p)%(2*p)-p
        #else:
            #perp_component = v - v.dot(coarseA[0])/(coarseA[0].dot(coarseA[0])) * coarseA[0]
            #if la.norm(perp_component) > tol:
                #coarseA[1] = v
                #if i < j:
                    #coarseb[1] = (b[i]+b[j]+p)%(2*p)-p
                #elif i == j:
                    #coarseb[1] = b[i]
                #elif i > j:
                    #coarseb[1] = (b[i]-b[abs(j)]+p)%(2*p)-p
                #break

    #if la.norm(coarseA[1]) == 0:
        #raise Exception("Poorly constructed A matrix: cannot find a pair of orthogonal vectors")

    ## find coarse solutions
    #try:
        #icA = la.inv(coarseA)
    #except:
        #raise Exception("Poorly constructed coarseA matrix: %s."%(coarseA))
    #try:
        #iA = la.inv(A.transpose().dot(A))
    #except:
        #raise Exception("Poorly constructed A matrix: %s."%(A.transpose().dot(A)))

    #if b.ndim > 2:
        #extra_shape = b.shape[1:]
        #flat_extra_dim = 1
        #for i in range(1, b.ndim):
            #flat_extra_dim = flat_extra_dim * b.shape[i]
        #coarseb.shape = (2, flat_extra_dim)
        #b.shape = (len(b), flat_extra_dim)
    #else:
        #extra_shape = None

    #result = icA.dot(coarseb)
    #for i in range(niter):
        #result = result + iA.dot(A.transpose().dot((b - A.dot(result) + p)%(2*p)-p))

    #if extra_shape is not None:
        #result.shape = tuple(np.concatenate(([2], extra_shape)))

    #return result

def extract_crosspol_ubl(data, info):
    '''input data should be xy/yx (2,...,bl)'''
    if len(data) != 2:
        raise AttributeError('Datas first demension need to have length 2 corresponding to xy/yx. Current input shape %s.'%data.shape)

    output_shape = np.array(data.shape)
    output_shape[-1] = info['nUBL']
    output = np.empty(output_shape, dtype='complex64')
    chisq = np.zeros(output_shape[:-1], dtype='float32')

    for u in range(info['nUBL']):
        blindex = info['subsetbl'][info['crossindex'][info['ublindex'][u][:,2].astype(int)]]
        ureversed = info['reversed'][info['ublindex'][u][:,2].astype(int)] == -1
        nreversed = np.sum(ureversed)
        if nreversed == 0:#no reversed
            output[..., u] = np.mean(data[..., blindex], axis=-1)
            chisq += np.linalg.norm(output[..., u][...,None] - data[..., blindex], axis=-1)**2
        elif nreversed == info['ublcount'][u]:
            output[..., u] = np.conjugate(np.mean(data[::-1, ..., blindex], axis=-1))
            chisq += np.linalg.norm(output[..., u][...,None] - np.conjugate(data[::-1, ..., blindex]), axis=-1)**2
        else:
            output[..., u] = (np.mean(data[..., blindex[~ureversed]], axis=-1) * (info['ublcount'][u] - nreversed) + np.conjugate(np.mean(data[::-1, ..., blindex[ureversed]], axis=-1)) * nreversed) / info['ublcount'][u]
            chisq += np.linalg.norm(output[..., u][...,None] - data[..., blindex[~ureversed]], axis=-1)**2 + np.linalg.norm(output[..., u][...,None] - np.conjugate(data[::-1, ..., blindex[ureversed]]), axis=-1)**2
    return output, chisq

# XXX data compression stuff belongs in another file
def deconvolve_spectra(spectra, window, band_limit, correction_weight=1e-15):
    '''solve for band_limit * 2 -1 bins, returns the deconvolved solution and 
    the norm of fitting error. All fft will be along first axis of spectra. 
    Input and outputs are in fourier space, window in real space'''
    if len(spectra) != len(window):
        raise ValueError("Input spectra and window function have unequal lengths %i %i."%(len(spectra), len(window)))
    #if np.sum(window) <= 2* band_limit - 1:
        #return np.zeros(2*band_limit - 1, dtype=np.array(spectra).dtype), np.inf
    fwindow = np.fft.fft(window) / len(window)
    band_limit_pass = np.zeros(len(fwindow), dtype='bool')
    band_limit_pass[:band_limit] = True
    if band_limit > 1:
        band_limit_pass[-(band_limit-1):] = True

    m = la.toeplitz(fwindow, np.roll(fwindow[::-1], 1)).astype('complex128')[:, band_limit_pass]
    mmi = la.inv(m.transpose().conjugate().dot(m) + np.identity(m.shape[1])*correction_weight)
    deconv_fdata = mmi.dot(m.transpose().conjugate()).dot(spectra)
    model_fdata = m.dot(deconv_fdata)
    return deconv_fdata, np.linalg.norm(model_fdata-spectra, axis = 0)

# XXX data compression stuff belongs in another file
def deconvolve_spectra2(spectra, window, band_limit, var=None, correction_weight=1e-15, correction_weight2=1e6):
    '''solve for band_limit * 2 -1 bins, returns the deconvolved solution 
    and the norm of fitting error. All fft will be along first axis of 
    spectra. Input and outputs are in real space, window also in real space'''
    if len(spectra) != len(window):
        raise ValueError("Input spectra and window function have unequal lengths %i %i."%(len(spectra), len(window)))
    #if np.sum(window) <= 2* band_limit - 1:
        #return np.zeros(2*band_limit - 1, dtype=np.array(spectra).dtype), np.inf
    if var is None:
        var = np.ones(len(window))
    elif len(var) != len(window):
        raise ValueError("Input var and window function have unequal lengths %i %i."%(len(var), len(window)))

    if np.sum(window) == 0:
        return np.zeros([2* band_limit - 1] + list(spectra.shape[1:])), np.zeros(spectra.shape[1:]), np.zeros_like(spectra), np.zeros((2* band_limit - 1, 2* band_limit - 1))+np.inf

    fwindow = np.fft.fft(window) / len(window)
    band_limit_pass = np.zeros(len(fwindow), dtype='bool')
    band_limit_pass[:band_limit] = True
    if band_limit > 1:
        band_limit_pass[-(band_limit-1):] = True

    #m = la.inv(la.dft(len(window))).dot(la.toeplitz(fwindow, np.roll(fwindow[::-1], 1)).astype('complex128')[:, band_limit_pass].dot(la.dft(2*band_limit - 1)))
    m = np.fft.ifft(la.toeplitz(fwindow, np.roll(fwindow[::-1], 1)).astype('complex128')[:, band_limit_pass].dot(la.dft(2*band_limit - 1)), axis=0)
    Ni = 1./np.copy(var)
    Ni[window==0] = np.nanmax(Ni) * correction_weight2

    #Ni = np.ones_like(window) + (1-window) * correction_weight2
    mmi = la.inv((m.transpose().conjugate() * Ni).dot(m) + np.identity(m.shape[1])*np.median(Ni[window==1])*correction_weight)
    mult_window = np.copy(window)
    mult_window.shape = [len(mult_window)] + [1]*(spectra.ndim-1)
    deconv_fdata = (mmi.dot(m.transpose().conjugate()) * Ni).dot(spectra*mult_window)
    model_fdata = m.dot(deconv_fdata)
    return deconv_fdata, np.linalg.norm(model_fdata-spectra*mult_window, axis = 0), model_fdata-spectra*mult_window, mmi

#  _____ _               
# |_   _(_)_ __  ___ _ _ 
#   | | | | '  \/ -_) '_|
#   |_| |_|_|_|_\___|_|  

# XXX utility stuff belongs in another file
class Timer:
    '''XXX DOCSTRING'''
    def __init__(self):
        self.time = time.time()
        self.start_time = self.time
        self.last_msg = None
        self.repeat_msg = 0

    def tick(self, msg='', mute=False):
        '''XXX DOCSTRING'''
        msg = str(msg)
        t = (float(time.time() - self.time)/60.)
        m = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
        if msg == self.last_msg:
            self.repeat_msg += 1
            if not mute:
                print msg + '*' + str(self.repeat_msg), "time elapsed: %f min"%t,
        else:
            self.repeat_msg = 0
            self.last_msg = msg
            if not mute:
             print msg, "Time elapsed: %f min."%t,
        if not mute:
            print "Memory usage 0: %.3fMB."%m
        sys.stdout.flush()
        self.time = time.time()
        return t, m

