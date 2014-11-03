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
    from scipy.stats import nanmedian

FILENAME = "calibration_omni.py"
julDelta = 2415020.# =julian date - pyephem's Observer date

infokeys = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
binaryinfokeys=['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B']


int_infokeys = ['nAntenna','nUBL','nBaseline']
intarray_infokeys = ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
float_infokeys = ['antloc','ubl','degenM','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
def read_redundantinfo_txt(infopath, verbose = False):
    METHODNAME = "read_redundantinfo_txt"
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


def write_redundantinfo_txt(info, infopath, overwrite = False, verbose = False):
    METHODNAME = "*write_redundantinfo_txt*"
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
        print FILENAME + "*" + METHODNAME + " MSG:", "Info file successfully written to %s. Time taken: %f minutes."%(infopath, (time.time()-timer)/60.)
    return


def write_redundantinfo(info, infopath, overwrite = False, verbose = False):
    METHODNAME = "*write_redundantinfo*"
    timer = time.time()
    if (not overwrite) and os.path.isfile(infopath):
        raise Exception("Error: a file exists at " + infopath + ". Use overwrite = True to overwrite.")
        return
    if (overwrite) and os.path.isfile(infopath):
        os.remove(infopath)
    marker = 9999999
    datachunk = [0 for i in range(len(binaryinfokeys)+1)]
    count = 0
    datachunk[count] = np.array([marker])         #start with a marker
    count += 1
    for key in binaryinfokeys:
        if key in ['antloc', 'ubl','degenM', 'AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']:  #'antloc',
            add = np.append(np.array(info[key]).flatten(),[marker])
            datachunk[count] = add
            count += 1
        elif key == 'ublindex':
            add = np.append(np.vstack(info[key]).flatten(),[marker])
            datachunk[count] = add
            count += 1
        elif key in ['A','B']:
            add = np.append(np.array(info[key].todense().flatten()).flatten(),[marker])
            datachunk[count] = add
            count += 1
        else:
            add = np.append(np.array(info[key]).flatten(),[marker])
            datachunk[count] = add
            count += 1
    datachunkarray = array('d',np.concatenate(tuple(datachunk)))
    outfile=open(infopath,'wb')
    datachunkarray.tofile(outfile)
    outfile.close()
    if verbose:
        print FILENAME + "*" + METHODNAME + " MSG:", "Info file successfully written to %s. Time taken: %f minutes."%(infopath, (time.time()-timer)/60.)
    return


def read_redundantinfo(infopath, verbose = False):
    METHODNAME = "read_redundantinfo"
    timer = time.time()
    if not os.path.isfile(infopath):
        raise Exception('Error: file path %s does not exist!'%infopath)
    with open(infopath) as f:
        farray=array('d')
        farray.fromstring(f.read())
        datachunk = np.array(farray)
        marker = 9999999
        markerindex=np.where(datachunk == marker)[0]
        rawinfo=np.array([np.array(datachunk[markerindex[i]+1:markerindex[i+1]]) for i in range(len(markerindex)-1)])

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

def importuvs(uvfilenames, totalVisibilityId, wantpols, nTotalAntenna = None, timingTolerance = math.pi/12/3600/100, init_mem = 4.e9, verbose = False):#tolerance of timing in radians in lst. init_mem is the initial memory it allocates for reading uv files.
    METHODNAME = "*importuvs*"
    ############################################################
    sun = ephem.Sun()
    #julDelta = 2415020
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
    deftime = int(init_mem / 8. / nfreq / (nant * (nant + 1) / 2))#use 4GB of memory by default.
    if verbose:
        print "Declaring initial array shape (%i, %i, %i, %i)..."%(deftime, len(wantpols), nant * (nant + 1) / 2, nfreq),
    sys.stdout.flush()
    try:
        data = np.zeros((deftime, len(wantpols), nant * (nant + 1) / 2, nfreq), dtype = 'complex64')
    except MemoryError:
        raise Exception("Failed to allocate %.2fGB of memory. Set init_mem keyword in Bytes for importuvs() to decrease initial memory allocation."%(init_mem/1.074e9))
    if verbose:
        print "Done."
    sys.stdout.flush()
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

def apply_omnigain_uvs(uvfilenames, omnigains, totalVisibilityId, info, wantpols, oppath, ano, adds=None, nTotalAntenna = None, overwrite = False, comment = '', verbose = False):
    METHODNAME = "*apply_omnigain_uvs*"
    ttotal = len(omnigains[wantpols.keys()[0]])
    ftotal = omnigains[wantpols.keys()[0]][0,0,3]
    if adds == None:
        adds = {}
        for key in wantpols.keys():
            adds[key] = np.zeros((ttotal, ftotal, len(totalVisibilityId)))
    if (ttotal != len(adds[wantpols.keys()[0]]) or ftotal != len(adds[wantpols.keys()[0]][0]) or len(totalVisibilityId) != len(adds[wantpols.keys()[0]][0,0])):
        raise Exception("Error: additives have different nTime or nFrequency or number of baseline!")
    if len(info) != len(omnigains) or len(info) != len(wantpols):
        raise Exception("Error: info and calparfilenames have different number of polarizations!")

    ####get some info from the first uvfile
    uv=ap.miriad.UV(uvfilenames[0])
    nfreq = uv.nchan;
    if nfreq != ftotal:
        raise Exception("Error: uv file %s and omnigains have different nFrequency!"%uvfilenames[0])
    if nTotalAntenna == None:
        nant = uv['nants'] # 'nants' should be the number of dual-pol antennas. PSA32 has a bug in double counting
    else:
        nant = nTotalAntenna

    if nant * (nant + 1) / 2 < len(totalVisibilityId):
        raise Exception("FATAL ERROR: Total number of antenna %d implies %d baselines whereas the length of totalVisibilityId is %d."%(nant, nant * (nant + 1) / 2, len(totalVisibilityId)))
    startfreq = uv['sfreq']
    dfreq = uv['sdf']
    del(uv)

    #######compute bl1dmatrix####each entry is 1 indexed with minus meaning conjugate, the bl here is the number in totalVisibilityId
    bl1dmatrix = np.zeros((nant, nant), dtype = 'int32')

    for a1a2, bl in zip(totalVisibilityId, range(len(totalVisibilityId))):
        a1, a2 = a1a2
        bl1dmatrix[a1, a2] = bl + 1
        bl1dmatrix[a2, a1] = - (bl + 1)
    ####load calpar from omnigain
    calpars = {}#bad antenna included
    for key in wantpols.keys():
        calpars[key] = (1. + np.zeros((ttotal, nant, nfreq),dtype='complex64'))
        calpars[key][:,info[key]['subsetant'],:] = omnigains[key][:,:,4::2] + 1.j * omnigains[key][:,:,5::2]


    #########start processing#######################
    t = []
    timing = []
    #datapulled = False
    for uvfile in uvfilenames:
        uvi = ap.miriad.UV(uvfile)
        if len(timing) > 0:
            if verbose:
                print FILENAME + METHODNAME + "MSG:", uvfile + ' after', timing[-1]#uv.nchan
                sys.stdout.flush()

        if oppath == None:
            oppath = os.path.abspath(os.path.dirname(os.path.dirname(uvfile + '/'))) + '/'
        opuvname = oppath + os.path.basename(os.path.dirname(uvfile+'/')) + ano + 'O'
        if verbose:
            print FILENAME + METHODNAME + "MSG: Creating %s"%opuvname
        if overwrite and os.path.isdir(opuvname):
            shutil.rmtree(opuvname)
        uvo = ap.miriad.UV(opuvname, status='new')
        uvo.init_from_uv(uvi)
        historystr = "Applied OMNICAL on %s: "%time.asctime(time.localtime(time.time()))
        uvo['history'] += historystr + comment + "\n"
        for preamble, data, flag in uvi.all(raw=True):
            uvo.copyvr(uvi)
            if len(t) < 1 or t[-1] != preamble[1]:#first bl of a timeslice
                t += [preamble[1]]

                if len(t) > ttotal:
                    raise Exception(FILENAME + METHODNAME + " MSG: FATAL ERROR: omnigain input array has length", omnigains[0].shape, "but the total length is exceeded when processing " + uvfile + " Aborted!")
            polwanted = False
            for pol in wantpols.keys():
                if wantpols[pol] == uvi['pol']:
                    a1, a2 = preamble[2]
                    bl = bl1dmatrix[a1, a2]
                    if bl > 0:
                        additive = adds[pol][len(t) - 1, :, bl - 1]
                    elif bl < 0:
                        additive = adds[pol][len(t) - 1, :, - bl - 1].conjugate()
                    else:
                        additive = 0
                        flag[:] = True
                    #print data.shape, additive.shape, calpars[pol][len(t) - 1, a1].shape
                    uvo.write(preamble, (data-additive)/calpars[pol][len(t) - 1, a1].conjugate()/calpars[pol][len(t) - 1, a2], flag)
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

def apply_omnical_uvs(uvfilenames, calparfilenames, totalVisibilityId, info, wantpols, oppath, ano, additivefilenames = None, nTotalAntenna = None, comment = '', overwrite= False):
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
        #for cpfn, adfn in zip(calparfilenames, additivefilenames):
            #historystr += os.path.abspath(cpfn) + ' ' + os.path.abspath(adfn) + ' '
        uvo['history'] += historystr + comment + "\n"
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
    def __init__(self, info, verbose=False):
        _O.RedundantInfo.__init__(self)
        if type(info) == type('a'):
            info = read_redundantinfo(info)
        elif type(info) != type({}):
            raise Exception("Error: info argument not recognized. It must be of either dictionary type (an info dictionary) *OR* string type (path to the info file).")
        if verbose:
            print "Converting info:",
            sys.stdout.flush()
        for key in info.keys():
            if verbose:
                print key,
                sys.stdout.flush()
            try:
                if key in ['At','Bt']:
                    tmp = []
                    nonzeros = np.array(info[key].nonzero()).transpose()
                    for i,j in nonzeros:
                        tmp += [[i, j, info[key][i,j]]]
                    #for i in range(info[key].shape[0]):
                        #for j in range(info[key].shape[1]):
                            #if info[key][i,j] != 0:
                                #tmp += [[i, j, info[key][i,j]]]
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
        if verbose:
            print "Done."
            sys.stdout.flush()

    def __getattribute__(self, key):
        try:
            if key in ['A','B']:
                #print key
                return sps.csr_matrix(_O.RedundantInfo.__getattribute__(self, key))
            elif key in ['At','Bt']:
                tmp = _O.RedundantInfo.__getattribute__(self, key+'sparse')
                matrix = np.zeros((self.nAntenna + self.nUBL, len(self.crossindex)))
                for i in tmp:
                    matrix[i[0],i[1]] = i[2]
                return sps.csr_matrix(matrix)
            elif key in ['ublindex']:
                ublindex = []
                for i in _O.RedundantInfo.__getattribute__(self, key):
                    while len(ublindex) < i[0] + 1:
                        ublindex.append(np.zeros((1,3)))
                    while len(ublindex[i[0]]) < i[1] + 1:
                        ublindex[i[0]] = np.array(ublindex[i[0]].tolist() + [[0,0,0]])
                    ublindex[i[0]][i[1]] = np.array(i[2:])
                return ublindex

            else:
                return _O.RedundantInfo.__getattribute__(self, key)
        except:
            raise Exception("Error retrieving %s item."%key)


    def get_info(self):
        info = {}
        for key in infokeys:
            try:
                #if key in ['A','B']:
                    ##print key
                    #info[key] = sps.csr_matrix(self.__getattribute__(key))
                #elif key in ['At','Bt']:
                    #tmp = self.__getattribute__(key+'sparse')
                    #matrix = np.zeros((info['nAntenna'] + info['nUBL'], len(info['crossindex'])))
                    #for i in tmp:
                        #matrix[i[0],i[1]] = i[2]
                    #info[key] = sps.csr_matrix(matrix)
                #elif key in ['ublindex']:
                    #ublindex = []
                    #for i in self.__getattribute__(key):
                        #while len(ublindex) < i[0] + 1:
                            #ublindex.append(np.zeros((1,3)))
                        #while len(ublindex[i[0]]) < i[1] + 1:
                            #ublindex[i[0]] = np.array(ublindex[i[0]].tolist() + [[0,0,0]])
                        #ublindex[i[0]][i[1]] = np.array(i[2:])
                    #info[key] = ublindex

                #else:
                    ##print key
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
        self.badUBLpair = []
        self.totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(self.nTotalAnt)])#PAPER miriad convention by default

        self.Info = None
        self.removeDegeneracy = True
        self.removeAdditive = False
        self.removeAdditivePeriod = -1
        self.convergePercent = 0.01 #convergence criterion in relative change of chi^2. By default it stops when reaches 0.01, namely 1% decrease in chi^2.
        self.maxIteration = 50 #max number of iterations in lincal
        self.stepSize = 0.3 #step size for lincal. (0, 1]. < 0.4 recommended.
        self.computeUBLFit = True

        self.nTime = 0
        self.nFrequency = 0

        self.utctime = None
        self.rawCalpar = None
        self.omnichisq = None
        self.omnigain = None
        self.omnifit = None

        if info != None:
            if type(info) == type({}):

                self.Info = RedundantInfo(info)
            elif type(info) == type('a'):
                self.read_redundantinfo(info)
            else:
                raise Exception(self.className + methodName + "Error: info argument not recognized. It must be of either dictionary type (an info dictionary) *OR* string type (path to the info file).")

    def read_redundantinfo(self, infopath, verbose = False):#redundantinfo is necessary for running redundant calibration. The text file should contain 29 lines each describes one item in the info.
        self.Info = RedundantInfo(read_redundantinfo(infopath, verbose = verbose), verbose = verbose)

    def write_redundantinfo(self, infoPath, overwrite = False, verbose = False):
        methodName = '.write_redundantinfo.'
        write_redundantinfo(self.Info.get_info(), infoPath, overwrite = overwrite, verbose = verbose)

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
        
        if len(np.array(rawinfo[1]).astype(int))%2 != 0 or min(np.array(rawinfo[1]).astype(int)) < 0:
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
                raise Exception(self.className + methodName + "Error: Format error in " + arrayinfopath + ": The antenna locations should start on the 4th line, with 3 numbers in each line!")
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
            print "Bad UBL indices:", self.badUBL


    def lincal(self, data, additivein, verbose = False):
        if data.ndim != 3 or data.shape[-1] != len(self.totalVisibilityId):
            raise Exception("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline(%i)"%len(self.totalVisibilityId))
        if data.shape != additivein.shape:
            raise Exception("Data shape error: data and additive in have different shapes.")
        self.nTime = len(data)
        self.nFrequency = len(data[0])
        if self.rawCalpar.shape != (len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)):
            raise Exception("ERROR: lincal called without a properly shaped self.rawCalpar! Excpeted shape is (%i, %i, %i)!"%(len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)))
        return _O.redcal(data[:,:,self.Info.subsetbl], self.rawCalpar, self.Info, additivein[:,:,self.Info.subsetbl], removedegen = int(self.removeDegeneracy), uselogcal = 0, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit))
        ##self.chisq = self.rawCalpar[:, :, 2]
        ##self.calpar = np.zeros((len(self.rawCalpar), len(self.rawCalpar[0]), self.nTotalAnt), dtype='complex64')
        ##self.calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
        ##self.bestfit = self.rawCalpar[:, :, (3 + 2 * self.Info.nAntenna):: 2] + 1.j * self.rawCalpar[:, :, (4 + 2 * self.Info.nAntenna):: 2]

    def logcal(self, data, additivein, verbose = False):
        if data.ndim != 3 or data.shape[-1] != len(self.totalVisibilityId):
            raise Exception("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline(%i)"%len(self.totalVisibilityId))
        if data.shape != additivein.shape:
            raise Exception("Data shape error: data and additive in have different shapes.")
        self.nTime = len(data)
        self.nFrequency = len(data[0])
        self.rawCalpar = np.zeros((len(data), len(data[0]), 3 + 2 * (self.Info.nAntenna + self.Info.nUBL)), dtype = 'float32')
        return _O.redcal(data[:,:,self.Info.subsetbl], self.rawCalpar, self.Info, additivein[:,:,self.Info.subsetbl], removedegen = int(self.removeDegeneracy), uselogcal = 1, maxiter=int(self.maxIteration), conv=float(self.convergePercent), stepsize=float(self.stepSize), computeUBLFit = int(self.computeUBLFit))

    def get_calibrated_data(self, data, additivein = None):
        if data.ndim != 3 or data.shape != (self.nTime, self.nFrequency, len(self.totalVisibilityId)):
            raise Exception("Data shape error: it must be a 3D numpy array of dimensions time * frequency * baseline (%i, %i, %i)"%(self.nTime, self.nFrequency, len(self.totalVisibilityId)))
        if additivein!= None and data.shape != additivein.shape:
            raise Exception("Data shape error: data and additivein have different shapes.")
        if data.shape[:2] != self.rawCalpar.shape[:2]:
            raise Exception("Data shape error: data and self.rawCalpar have different first two dimensions.")

        calpar = np.ones((len(self.rawCalpar), len(self.rawCalpar[0]), self.nTotalAnt), dtype='complex64')
        calpar[:,:,self.Info.subsetant] = (10**(self.rawCalpar[:, :, 3: (3 + self.Info.nAntenna)])) * np.exp(1.j * self.rawCalpar[:, :, (3 + self.Info.nAntenna): (3 + 2 * self.Info.nAntenna)])
        if additivein == None:
            return apply_calpar(data, calpar, self.totalVisibilityId)
        else:
            return apply_calpar(data - additivein, calpar, self.totalVisibilityId)


    def get_omnichisq(self):
        if self.utctimes == None or self.rawCalpar == None:
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
        if self.utctimes == None or self.rawCalpar == None:
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
        if self.utctimes == None or self.rawCalpar == None:
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



    def diagnose(self, data = None, additiveout = None, verbose = True, healthbar = 2, ubl_healthbar = 50, warn_low_redun = False):
        errstate = np.geterr()
        np.seterr(invalid = 'ignore')
        checks = 1
        bad_count = np.zeros((3,self.Info.nAntenna), dtype='int')
        bad_ubl_count = np.zeros(self.Info.nUBL, dtype='int')
        median_level = nanmedian(nanmedian(self.rawCalpar[:,:,3:3+self.Info.nAntenna], axis= 0), axis= 1)
        bad_count[0] = np.array([(np.abs(self.rawCalpar[:,:,3+a] - median_level) >= .5).sum() for a in range(self.Info.nAntenna)])**2

        if data != None and data.shape[:2] == self.rawCalpar.shape[:2]:
            checks += 1
            subsetbl = self.Info.subsetbl
            crossindex = self.Info.crossindex
            ncross = len(self.Info.crossindex)
            bl1dmatrix = self.Info.bl1dmatrix
            ant_level = np.array([np.median(np.abs(data[:, :, [subsetbl[crossindex[bl]] for bl in bl1dmatrix[a] if bl < ncross]]), axis = 2) for a in range(self.Info.nAntenna)])
            median_level = np.median(ant_level, axis = 0)
            bad_count[1] = np.array([(np.abs(ant_level[a] - median_level)/median_level >= .667).sum() for a in range(self.Info.nAntenna)])**2

        if additiveout != None and additiveout.shape[:2] == self.rawCalpar.shape[:2]:
            checks += 1
            subsetbl = self.Info.subsetbl
            crossindex = self.Info.crossindex
            ncross = len(self.Info.crossindex)
            bl1dmatrix = self.Info.bl1dmatrix
            ant_level = np.array([np.median(np.abs(additiveout[:, :, [crossindex[bl] for bl in bl1dmatrix[a] if bl < ncross]]), axis = 2) for a in range(self.Info.nAntenna)])
            median_level = np.median(ant_level, axis = 0)
            bad_count[2] = np.array([(np.abs(ant_level[a] - median_level)/median_level >= .667).sum() for a in range(self.Info.nAntenna)])**2

            ublindex = [np.array(index).astype('int')[:,2] for index in self.Info.ublindex]
            ubl_level = np.array([np.median(np.abs(additiveout[:, :, [crossindex[bl] for bl in ublindex[u]]]), axis = 2) for u in range(self.Info.nUBL)])
            median_level = np.median(ubl_level, axis = 0)
            bad_ubl_count += np.array([((ubl_level[u] - median_level)/median_level >= .667).sum() for u in range(self.Info.nUBL)])**2
            #print median_level
        np.seterr(invalid = errstate['invalid'])
        bad_count = (np.mean(bad_count,axis=0)/float(self.nTime * self.nFrequency)**2 * 100).astype('int')
        bad_ubl_count = (bad_ubl_count/float(self.nTime * self.nFrequency)**2 * 100).astype('int')
        if verbose:
            #print bad_ant_cnt, bad_ubl_cnt
            print "DETECTED BAD ANTENNA ABOVE HEALTH THRESHOLD %i: "%healthbar
            for a in range(len(bad_count)):
                if bad_count[a] > healthbar:
                    print "antenna #%i, vector = %s, badness = %i"%(self.Info.subsetant[a], self.Info.antloc[a], bad_count[a])
            #print ""
            if additiveout != None and additiveout.shape[:2] == self.rawCalpar.shape[:2] and ubl_healthbar != 100:
                print "DETECTED BAD BASELINE TYPE ABOVE HEALTH THRESHOLD %i: "%ubl_healthbar
                for a in range(len(bad_ubl_count)):
                    if bad_ubl_count[a] > ubl_healthbar and (self.Info.ublcount[a] > 5 or (warn_low_redun)):
                        print "index #%i, vector = %s, redundancy = %i, badness = %i"%(a, self.Info.ubl[a], self.Info.ublcount[a], bad_ubl_count[a])
                #print ""
        return bad_count, bad_ubl_count


    def compute_redundantinfo(self, arrayinfoPath = None):
        if arrayinfoPath != None and os.path.isfile(arrayinfoPath):
            self.read_arrayinfo(arrayinfoPath)
        if np.linalg.norm(self.antennaLocation) == 0:
            raise Exception("Error: compute_redundantinfo() called before self.antennaLocation is specified. Use configFilePath option when calling compute_redundantinfo() to specify array info file, or manually set self.antennaLocation for the RedundantCalibrator instance.")
        #timer = Timer()
        #nAntenna and subsetant : get rid of the bad antennas
        nant=len(self.antennaLocation)
        subsetant=[i for i in range(nant) if i not in self.badAntenna]
        nAntenna=len(subsetant)
        antloc=[self.antennaLocation[ant] for ant in subsetant]
        #timer.tick('a')
        ##########################################################################################
        #find out ubl
        #use the function compute_UBL to find the ubl
        tolerance=self.antennaLocationTolerance;
        ublall=self.compute_UBL(tolerance)
        #timer.tick('b')
        #################################################################################################
        #calculate the norm of the difference of two vectors (just la.norm actually)
        def dis(a1,a2):
            return np.linalg.norm(np.array(a1)-np.array(a2))
        #find badUBL with badUBLpair
        def find_ublindex_all(pair):
            for i in range(len(ublall)):
                if dis(self.antennaLocation[pair[0]]-self.antennaLocation[pair[1]],ublall[i]) < tolerance or dis(self.antennaLocation[pair[0]]-self.antennaLocation[pair[1]],-ublall[i]) < tolerance:
                    return i
            return None
            #raise Exception("Error: something wrong in identifying badUBL from badUBLpair")    #delete this line afterwards
        for p in self.badUBLpair:
            self.badUBL.append(find_ublindex_all(p))
        self.badUBL = [i for i in self.badUBL if i != None]
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

        #correct the orders of pairs in goodpair
        def correct_pairorder(pair):
            try:
                self.totalVisibilityId.tolist().index([pair[0],pair[1]])
                return True
            except:
                try:
                    self.totalVisibilityId.tolist().index([pair[1], pair[0]])
                    return False
                except:
                    return None
                    
        #exclude pairs that are not in totalVisibilityId
        temp = []
        for p in goodpairs:
            cond = correct_pairorder([subsetant[p[0]],subsetant[p[1]]])
            if cond == True:
                temp.append(p)
            if cond == False:
                temp.append(p[::-1])
        goodpairs = temp
        
        #goodpairs = [correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) for p in goodpairs if (correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) != None and correct_pairorder([subsetant[p[0]],subsetant[p[1]]]) == True)]  #correct_pairorder([subsetant[p[0]],subsetant[p[1]]])
        nBaseline=len(goodpairs)

        #from a pair of good antenna index to baseline index
        subsetbl = np.array([self.get_baseline([subsetant[bl[0]],subsetant[bl[1]]]) for bl in goodpairs])
        #timer.tick('c')
        ##################################################################################
        #bltoubl: cross bl number to ubl index
        def findublindex(pair,ubl=ubl):
            i=pair[0]
            j=pair[1]
            for k in range(len(ubl)):
                if dis(antloc[i]-antloc[j],ubl[k])<tolerance or dis(antloc[i]-antloc[j],-ubl[k])<tolerance:
                    return k
            print pair
            return "no match"
        bltoubl=[];
        for p in goodpairs:
            if p[0]!=p[1]:
                bltoubl.append(findublindex(p))
        #timer.tick('d')
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
                print "something's wrong with bltoubl"
        #timer.tick('e')
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
        #timer.tick('f')
        ###################################################
        #bl2d:  from 1d bl index to a pair of antenna numbers
        bl2d=[]
        for pair in goodpairs:
            bl2d.append(pair)#(pair[::-1])
        bl2d=np.array(bl2d)
        #timer.tick('g')
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
        #timer.tick('h')
        ####################################################################################
        #ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
        countdict={}
        for bl in bltoubl:
            countdict[bl]=[]

        for i in range(len(crosspair)):
            ant1=crosspair[i][0]
            ant2=crosspair[i][1]
            countdict[bltoubl[i]].append([ant1,ant2,i])  #([ant1,ant2,i])

        ublindex=[]
        for i in range(nUBL):
            ublindex.append(countdict[i])
        #turn each list in ublindex into np array
        for i in range(len(ublindex)):
            ublindex[i]=np.array(ublindex[i])
        ublindex=np.array(ublindex)
        #timer.tick('i')
        ###############################################################################
        #bl1dmatrix: a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
                #I suppose 99999 for bad and auto baselines?
        bl1dmatrix=99999*np.ones([nAntenna,nAntenna],dtype='int16')
        for i in range(len(crosspair)):
            bl1dmatrix[crosspair[i][1]][crosspair[i][0]]=i
            bl1dmatrix[crosspair[i][0]][crosspair[i][1]]=i
        #timer.tick('j')
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
            B[i][crosspair[i][0]]=reverse[i]*-1   #1
            B[i][crosspair[i][1]]=reverse[i]*1  #-1
            B[i][nAntenna+bltoubl[i]]=1
        B=sps.csr_matrix(B)
        #timer.tick('k')
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
        #timer.tick('l')
        self.Info = RedundantInfo(info)
        #timer.tick('m')


    #inverse function of totalVisibilityId, calculate the baseline index from the antenna pair. It allows flipping of a1 and a2, will return same result
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
        try:
            return self.totalVisibilityId.tolist().index([pair[0],pair[1]])
        except:
            try:
                return self.totalVisibilityId.tolist().index([pair[1], pair[0]])
            except:
                #raise Exception("Error: antenna pair %s not found in self.totalVisibilityId."%pair)
                return None
        #Eric's old code. It's buggy and assumes totalVisibilityId contains a1,a2 where a2<=a1 always
        #sortp = np.array(sorted(pair))
        #for i in range(len(self.totalVisibilityId)):
            #if self.totalVisibilityId[i][0] == sortp[0] and self.totalVisibilityId[i][1] == sortp[1]:
                #return i
        #raise Exception("antenna index out of range")

                
                
    #compute_UBL returns the average of all baselines in that ubl group
    def compute_UBL_old(self,tolerance = 0.1):
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
        
        
    #compute_UBL returns the average of all baselines in that ubl group
    def compute_UBL(self,tolerance = 0.1):
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





def omniview(data, info, plotrange = None, title = ''):
    import matplotlib.pyplot as plt
    d=data[info['subsetbl']][info['crossindex']]
    if plotrange == None:
        plotrange = 1.2*np.nanmax(np.abs(d))
    ubl = 0
    colors=[]
    colorgrid = int(math.ceil((info['nUBL']/12.+1)**.34))
    for red in range(colorgrid):
        for green in range(colorgrid):
            for blue in range(colorgrid):
                #print red, green, blue
                colors += [(np.array([red, green, blue]).astype('float')/(colorgrid - 1)).tolist()]
    #colors.remove([0,0,0])
    colors.remove([1,1,1])

    for marker in ["o", "v", "^", "<", ">", "8", "s", "p", "h", (6,1,0), (8,1,0), "d"]:
        for color in colors:
            #print info['ublindex'][ubl][:,2]
            #print marker, color
            plt.scatter(np.real(d[np.array(info['ublindex'][ubl][:,2]).astype('int')]),np.imag(d[np.array(info['ublindex'][ubl][:,2]).astype('int')])*info['reversed'][np.array(info['ublindex'][ubl][:,2]).astype('int')], marker=marker, color=color)
            ubl += 1
            if ubl == info['nUBL']:
                plt.xlabel('Real')
                plt.ylabel('Imag')
                plt.title(title)
                plt.grid(True)
                plt.axis([-plotrange, plotrange, -plotrange, plotrange])
                plt.axes().set_aspect('equal')
                plt.axes().text(-0.9*plotrange, -0.9*plotrange, "#Ant:%i\n#UBL:%i"%(info['nAntenna'],info['nUBL']),bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))
                plt.show()
                return

def lin_depend(v1, v2, tol = 0):#whether v1 and v2 are linearly dependent
    if len(v1) != len(v2):
        raise Exception("Length mismatch %i vs %i."%(len(v1), len(v2)))
    if la.norm(v1) == 0:
        return True
    return la.norm(np.dot(v1, v2)/np.dot(v1, v1) * v1 - v2) <= tol

def find_solution_path(info, rawcal_ubl=[], tol = 0.0, verbose=False):#return (intialantenna, solution_path) for raw calibration. solution path contains a list of [(a0, a1, crossubl), a] = [(ublindex entry), (which ant is solved, 0 or 1)]. When raw calibrating, initialize calpar to have [0] at initial antenna, then simply iterate along the solution_path, use crossubl and a0 or a1 specified by a to solve for the other a1 or a0 and append it to calpar. Afterwards, use mean angle on calpars
    ###select 2 ubl for calibration
    if rawcal_ubl == []:
        ublcnt_tmp = info['ublcount'].astype('float')
        rawcal_ubl += [np.argmax(ublcnt_tmp)]
        ublcnt_tmp[rawcal_ubl[-1]] = np.nan
        rawcal_ubl += [np.nanargmax(ublcnt_tmp)]
        ublcnt_tmp[rawcal_ubl[-1]] = np.nan
        #while np.allclose(info['ubl'][rawcal_ubl[0]]/(la.norm(info['ubl'][rawcal_ubl[0]])/la.norm(info['ubl'][rawcal_ubl[1]])), info['ubl'][rawcal_ubl[1]]) or np.allclose(info['ubl'][rawcal_ubl[0]]/(la.norm(info['ubl'][rawcal_ubl[0]])/la.norm(info['ubl'][rawcal_ubl[1]])), -info['ubl'][rawcal_ubl[1]]):
        while lin_depend(info['ubl'][rawcal_ubl[0]], info['ubl'][rawcal_ubl[1]], tol=tol):
            try:
                rawcal_ubl[1] = np.nanargmax(ublcnt_tmp)
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
    A = np.array([bl1[:2], bl2[:2]])
    remove_Matrix = (np.array(info['antloc'])- info['antloc'][initialant])[:,:2].dot(la.pinv(A.transpose().dot(A)).dot(A.transpose()))
    degeneracy_remove = [a1, a2, remove_Matrix]
    if verbose:
        print "Degeneracy: a1 = %i, a2 = %i"%(info['subsetant'][a1], info['subsetant'][a2])
    return initialant, solution_path, additional_solution_path, degeneracy_remove, (unsolved_ant == [])

def meanAngle(a):
    return np.angle(np.mean(np.exp(1.j*np.array(a))))
def medianAngle(a):
    return np.angle(np.median(np.cos(np.array(a))) + 1.j * np.median(np.sin(np.array(a))))

def raw_calibrate(data, info, initant, solution_path, additional_solution_path, degeneracy_remove):
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

class Timer():
    def __init__(self):
        self.time = time.time()

    def tick(self, msg=''):
        print msg, "time elapsed: %f min"%(float(time.time() - self.time)/60.)
        sys.stdout.flush()
        self.time = time.time()
        
def remove_one_antenna(Info,badant):
    info = Info.get_info()
    #nAntenna and antloc
    nAntenna = info['nAntenna']-1 
    badindex = list(info['subsetant']).index(badant)     #the index of the bad antenna in previous subsetant

    subsetant = list(info['subsetant'])[:]
    subsetant.pop(badindex)      #delete the bad antenna from subsetant
    antloc = np.delete(np.array(info['antloc']),badindex,0)   #delete the bad antenna from antloc
         
    #ubl and nUBL
    index = 0              #to keep track of the index of ubl the loop is at
    deletelist = []
    for ubl in info['ublindex']:    
        if len(ubl) > 1:
            index += 1
        elif ubl[0,0] == subsetant[badant] or ubl[0,1] == subsetant[badant] :
            deletelist.append(index)
            index += 1

    ubl = info['ubl'][:]
    ubl = np.array([ubl[i] for i in range(len(ubl)) if i not in deletelist])
    nUBL=len(ubl);

    #subsetbl and nBaseline     
    goodpairs_old = [i[::-1] for i in info['bl2d']]    #the old goodpairs
    goodpairs_index = [i for i in range(len(goodpairs_old)) if badindex not in goodpairs_old[i]]       #the index of goodpairs that doesn't contain the bad antenna
    temp = np.array([goodpairs_old[i] for i in goodpairs_index])   #the new goodpairs with antenna number (all antenna)
    goodpairs = np.zeros(temp.shape)
    for i in range(len(temp)):
        for j in range(len(temp[i])):
            if temp[i,j] > badindex:
                goodpairs[i,j] = temp[i,j]-1
            else:
                goodpairs[i,j] = temp[i,j]

    subsetbl = [info['subsetbl'][i] for i in goodpairs_index]  #the new subsetbl
    nBaseline = len(subsetbl)

    counter = 0
    ubl_old2new = np.zeros([len(info['ubl'])],dtype = 'int')     #from old ubl index to new ubl index
    for i in range(len(info['ubl'])):
        if i in deletelist:
            ubl_old2new[i] = counter
        else:
            ubl_old2new[i] = counter
            counter += 1

    bltoubl = []
    for i in range(len(info['crossindex'])):
        pair = [info['subsetant'][index] for index in info['bl2d'][info['crossindex'][i]]]   #get the pair of antenna from each crossindex 
        if badant in pair:
            pass
        else:
            bltoubl.append(ubl_old2new[info['bltoubl'][i]])   #append the new ubl index that doesn't have the bad antenna
    bltoubl = np.array(bltoubl)
    #################################################################################
    #reversed:   cross only bl if reversed -1, otherwise 1
    def dis(a1,a2):    #calculate the norm of the difference of two vectors
        return np.linalg.norm(np.array(a1)-np.array(a2))

    crosspair_old = []
    for p in goodpairs_old:
        if p[0]!=p[1]:
            crosspair_old.append(p)
    goodcross = []
    for i in range(len(crosspair_old)):
        if badindex not in crosspair_old[i]:
            goodcross.append(i)

    crosspair=[]
    for p in goodpairs:
        if p[0]!=p[1]:
            crosspair.append(p)
            
    reverse=[info['reversed'][i] for i in goodcross]
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
    ############################################################################
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
    return RedundantInfo(info)
