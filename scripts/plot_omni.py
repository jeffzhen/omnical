#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem, optparse, sys
import omnical.calibration_omni as omni
import cPickle as pickle
import scipy.signal as ss
import scipy.ndimage.filters as sfil
import matplotlib.pyplot as plt
FILENAME = "omnical_PSA128.py"
print "#Omnical Version %s#"%omni.__version__


######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

ap.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-t', '--time', action = 'store', type = 'int', default = 0, help = "Time for uv file plotting. Need to match uv['time'] entry.")
o.add_option('-f', '--frequency', action = 'store', type = 'int', default = 0, help = 'Frequency channel number for the uv file.')
o.add_option('-r', '--range', action = 'store', default = None, help = 'Plot range symmetric for single number or use min_max.')
o.add_option('-a', '--antenna', action = 'store', default = '0', help = 'antenna number or antenna pair for omnigain/omnifit.')

o.add_option('-i', '--infopath', action = 'store', default = None, help = 'Manually feed redundantinfo file to read. Only necessary for uv file plotting.')
o.add_option('-m', '--mode', action = 'store', default = 'phs', help = 'Specify plot mode for omnigain or omnifit: amp, phs, real, imag')



opts, args = o.parse_args(sys.argv[1:])
datafiles = [os.path.expanduser(arg) for arg in args]
if opts.infopath is not None:
    infopath = os.path.expanduser(opts.infopath)
else:
    infopath = None
time = int(opts.time)
frequency = int(opts.frequency)
antenna = opts.antenna
mode = opts.mode
plotrange = opts.range
if plotrange is not None:
    try:
        ranges = [float(v) for v in plotrange.split('_')]
        if len(ranges) == 1:
            assert ranges[0] > 0
            ranges = [-ranges[0], ranges[0]]
        else:
            assert len(ranges) == 2 and ranges[0] < ranges[1]
    except:
        raise IOError("Invalid range option %s."%opts.range)
    force_range = True
else:
    force_range = False

if '.omnigain' in datafiles[0]:
    plottype = 'omnigain'
elif '.omnifit' in datafiles[0]:
    plottype = 'omnifit'
elif '.omnichisq' in datafiles[0]:
    plottype = 'omnichisq'
elif '.uv' in datafiles[0] and os.path.isfile(datafiles[0] + '/visdata'):
    plottype = 'uv'
else:
    raise IOError("File format not recognized. Must contain .omnigain, .omnifit, .omnichisq, or .uv")

######################################################################
############## read redundant info and antenna numbers; perform sanity checks ###################################
######################################################################
if plottype == 'uv':
    if len(datafiles) != 1:
        raise IOError("For uv file plotting only 1 file accepted. Received %i."%len(datafiles))
    if not os.path.isdir(datafiles[0]):
        raise IOError("UV file not found: %s."%(datafiles[0]))
    if infopath is None:
        raise IOError("For uv file plotting you must supply a redundant info file.")
    elif not os.path.isfile(infopath):
        raise IOError("Redundant info file not found: %s."%infopath)
    else:
        info = omni.read_redundantinfo(infopath)

else:
    if mode not in ['amp', 'phs', 'real', 'imag']:
        raise IOError("Specified mode %s not recognized. Only allow amp, phs, real, imag.")

    for datafile in datafiles:
        if not os.path.isfile(datafile):
            raise IOError("%s file not found: %s."%(plottype, datafile))

    if infopath is None:
        infopath = datafiles[0].replace(plottype, 'binfo')
    if not os.path.isfile(infopath):
        raise IOError("Redundant info file not found: %s."%infopath)
    else:
        info = omni.read_redundantinfo(infopath)

    if plottype == 'omnigain':
        if '_' in antenna:
            raise IOError("Antenna option %s is not valid for %s file. Accepts single antenna or antenna numbers seperated by ','."%(antenna, plottype))
        else:
            try:
                antennas = np.array([int(a) for a in antenna.split(',')])
            except ValueError:
                raise IOError("Antenna option %s is not valid for %s file. Accepts single antenna or antenna numbers seperated by ','."%(antenna, plottype))

            try:
                antennas = np.array([list(info['subsetant']).index(a) for a in antennas])#good antenna index rather than the literall antenna index
            except ValueError:
                raise IOError("One of the antennas are not valid. The good antennas are: %s."%info['subsetant'])

    if plottype == 'omnifit':
        if ('_' not in antenna) or (antenna.count('_') != (antenna.count(',') + 1)):
            raise IOError("Antenna option %s is not valid for %s file. Accepts single antenna pair like 0_1 or many antenna pairs seperated by ','."%(antenna, plottype))
        else:
            antennas = np.array([[int(a) for a in apair.split('_')] for apair in antenna.split(',')])

            for apair in antennas:
                if len(apair) != 2:
                    raise IOError("Antenna option %s is not valid for %s file. Accepts single antenna pair like 0_1 or many antenna pairs seperated by ','. Got a bad pair of length %i."%(antenna, plottype, len(apair)))
            try:
                antennas = np.array([[list(info['subsetant']).index(a) for a in apair] for apair in antennas])#good antenna index rather than the literall antenna index
            except ValueError:
                raise IOError("One of the antennas are not valid. The good antennas are: %s."%info['subsetant'])

            ubls = np.zeros(len(antennas), dtype = int)
            for i, apair in enumerate(antennas):
                crossbl = info['bl1dmatrix'][apair[0]][apair[1]]
                if crossbl > len(info['bltoubl']) or crossbl < 0:
                    raise IOError("One of the antenna pairs %s is not valid."%apair)
                ubls[i] = info['bltoubl'][crossbl]


######################################################################
############## omnigain ###################################
######################################################################
if plottype == 'omnigain':
    p = 0
    nprefix = 2
    for datafile in datafiles:
        data = np.fromfile(datafile, dtype='complex64')
        nf = int(np.imag(data[nprefix - 1]))
        nt = len(data) / len(info['subsetant']) / (nf + nprefix)
        if len(data) != nt * len(info['subsetant']) * (nf + nprefix):
            raise IOError("File %s is incompatible with shape %s inferred from info file %s."%(datafile, (nt, len(info['subsetant']), (nf + nprefix)), infopath))
        else:
            data.shape = (nt, len(info['subsetant']), (nf + nprefix))

            if mode == 'amp':
                data = np.abs(data[:, antennas, nprefix:])
            elif mode == 'phs':
                data = np.angle(data[:, antennas, nprefix:])
            elif mode == 'real':
                data = np.real(data[:, antennas, nprefix:])
            elif mode == 'imag':
                data = np.imag(data[:, antennas, nprefix:])
        for i,a in enumerate(antennas):
            p = p + 1
            plt.subplot(len(antennas), len(datafiles), p)

            plotdata = data[:, i]
            if not force_range:
                ranges = [np.percentile(plotdata, 10), np.percentile(plotdata, 90)]
            plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
            plt.title('%s Ant#%i'%(os.path.basename(datafile), info['subsetant'][a]))
            plt.colorbar()
    plt.show()


######################################################################
############## omnifit ###################################
######################################################################
if plottype == 'omnifit':
    p = 0
    nprefix = 3
    for datafile in datafiles:
        data = np.fromfile(datafile, dtype='complex64')
        nf = int(np.imag(data[nprefix - 1]))
        nt = len(data) / info['nUBL'] / (nf + nprefix)
        if len(data) != nt * info['nUBL'] * (nf + nprefix):
            raise IOError("File %s is incompatible with shape %s inferred from info file %s."%(datafile, (nt, info['nUBL'], (nf + nprefix)), infopath))
        else:
            data.shape = (nt, info['nUBL'], (nf + nprefix))

            if mode == 'amp':
                data = np.abs(data[:, ubls, nprefix:])
            elif mode == 'phs':
                data = np.angle(data[:, ubls, nprefix:])
            elif mode == 'real':
                data = np.real(data[:, ubls, nprefix:])
            elif mode == 'imag':
                data = np.imag(data[:, ubls, nprefix:])
        for i,u in enumerate(ubls):
            p = p + 1
            plt.subplot(len(ubls), len(datafiles), p)

            plotdata = data[:, i]
            if not force_range:
                ranges = [np.percentile(plotdata, 5), np.percentile(plotdata, 95)]
            plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
            plt.title('%s, baseline [%.1f, %.1f, %.1f]'%(os.path.basename(datafile), info['ubl'][u][0], info['ubl'][u][1], info['ubl'][u][2]))
            plt.colorbar()
    plt.show()


######################################################################
############## omnichisq ###################################
######################################################################
if plottype == 'omnichisq':
    p = 0
    nprefix = 3
    for datafile in datafiles:
        data = np.fromfile(datafile, dtype='float32')
        nf = int(data[nprefix - 1])
        nt = len(data) / (nf + nprefix)
        if len(data) != nt * (nf + nprefix):
            raise IOError("File %s is not a valid omnichisq file."%datafile)
        else:
            data.shape = (nt, (nf + nprefix))
            plotdata = data[:, nprefix+1:]
        p = p + 1
        plt.subplot(1, len(datafiles), p)
        if not force_range:
            ranges = [np.percentile(plotdata, 5), np.percentile(plotdata, 95)]
        plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
        plt.title('%s, chi^2'%os.path.basename(datafile))
        plt.colorbar()
    plt.show()
######################################################################
############## omniview uv ###################################
######################################################################


