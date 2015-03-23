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

o.add_option('-t', '--time', action = 'store', type = 'int', default = None, help = "Time index for uv file plotting.")
o.add_option('-f', '--frequency', action = 'store', type = 'int', default = None, help = 'Frequency channel number for the uv file.')
o.add_option('-p', '--pol', action = 'store', default = 'xx', help = 'Polarization string for uv file plotting.')
o.add_option('-o', '--oppath', action = 'store', default = None, help = 'output path including file name.')
o.add_option('-r', '--range', action = 'store', default = None, help = 'Plot range symmetric for single number or use min_max.')
o.add_option('-a', '--antenna', action = 'store', default = '0', help = 'antenna number or antenna pair for omnigain/omnifit.')

o.add_option('-i', '--infopath', action = 'store', default = None, help = 'Manually feed redundantinfo file to read. Only necessary for uv file plotting.')
o.add_option('-m', '--mode', action = 'store', default = 'phs', help = 'Specify plot mode for omnigain or omnifit: amp, phs, real, imag')
o.add_option('--overwrite', action = 'store_true', help = 'whether to overwrite output file.')
o.add_option('-s', '--suppress', action = 'store_true', help = 'whether to suppress pop-out plot for scripting purpose.')
o.add_option('-u', '--unflag', action = 'store_true', help = 'whether to suppress flagging.')



opts, args = o.parse_args(sys.argv[1:])
if len(args) == 0:
    raise IOError("No files to plot.")
else:
    datafiles = [os.path.expanduser(arg) for arg in args]
if opts.infopath is not None:
    infopath = os.path.expanduser(opts.infopath)
else:
    infopath = None

overwrite = opts.overwrite
suppress = opts.suppress
if opts.oppath is not None:
    oppath = os.path.expanduser(opts.oppath)
    if os.path.isdir(oppath) or (os.path.dirname(oppath) != '' and not os.path.isdir(os.path.dirname(oppath))):
        raise IOError("%s is not a valid path for outputting the plot."%opts.oppath)
    if os.path.isfile(oppath) and not overwrite:
        raise IOError("%s exists. Use --overwrite to overwrite."%oppath)
else:
    oppath = None
    if suppress:
        raise IOError("You have asked to suppress pop-out plotting with -s/--supress without supplying a valid path to store the plot with -o.")

if opts.time is not None:
    time = int(opts.time)
else:
    time = None
if opts.frequency is not None:
    frequency = int(opts.frequency)
else:
    frequency = None

pol = opts.pol
antenna = opts.antenna
mode = opts.mode
unflag = opts.unflag
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

    if time is None or frequency is None:
        raise IOError("For uv file plotting you must supply both a time index (-t) and a frequency index (-f).")

else:
    if time is not None and frequency is not None:
        raise IOError("Cannot accept both -t and -f options for %s file. Accepts at most one out of the two options."%(plottype))

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
            reverse = np.zeros(len(antennas), dtype = bool)
            for i, apair in enumerate(antennas):
                crossbl = info['bl1dmatrix'][apair[0]][apair[1]]
                if crossbl > len(info['bltoubl']) or crossbl < 0:
                    raise IOError("One of the antenna pairs %s is not valid."%apair)
                ubls[i] = info['bltoubl'][crossbl]
                reverse[i] = info['reversed'][crossbl]


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
        data.shape = (nt, len(info['subsetant']), (nf + nprefix))

        if mode == 'amp':
            data = np.abs(data[:, antennas, nprefix:])
            if force_range:
                ranges[0] = max(ranges[0], 0)
        elif mode == 'phs':
            data = np.angle(data[:, antennas, nprefix:])
        elif mode == 'real':
            data = np.real(data[:, antennas, nprefix:])
        elif mode == 'imag':
            data = np.imag(data[:, antennas, nprefix:])

        flagpath = datafile.replace(plottype, 'omniflag')
        if os.path.isfile(flagpath) and not unflag:
            flag = np.fromfile(flagpath, dtype='bool').reshape((nt, nf))
        else:
            flag = np.zeros((nt, nf), dtype='bool')

        if time is not None:
            flag = flag[time]
        elif frequency is not None:
            flag = flag[:,frequency]

        for i,a in enumerate(antennas):
            p = p + 1
            plt.subplot(len(antennas), len(datafiles), p)
            if time is not None:
                plotdata = data[time, i]
            elif frequency is not None:
                plotdata = data[:, i, frequency]
            else:
                plotdata = data[:, i]

            if not force_range:
                ranges = [np.percentile(plotdata, 10), np.percentile(plotdata, 90)]

            plotdata[flag] = np.nan
            if len(plotdata.shape) == 2:
                plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
                plt.colorbar()
            else:
                plt.plot(plotdata)
                plt.ylim(ranges)
            plt.title('%s Ant#%i %s'%(os.path.basename(datafile), info['subsetant'][a], mode))

    if oppath is not None:
        plt.savefig(oppath, bbox_inches='tight')
    if not suppress:
        plt.show()
    else:
        plt.close()


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
        data.shape = (nt, info['nUBL'], (nf + nprefix))

        if mode == 'amp':
            data = np.abs(data[:, ubls, nprefix:])
            if force_range:
                ranges[0] = max(ranges[0], 0)
        elif mode == 'phs':
            data = np.angle(data[:, ubls, nprefix:])
            data[:,reverse] = -data[:,reverse]
        elif mode == 'real':
            data = np.real(data[:, ubls, nprefix:])
        elif mode == 'imag':
            data = np.imag(data[:, ubls, nprefix:])
            data[:,reverse] = -data[:,reverse]

        flagpath = datafile.replace(plottype, 'omniflag')
        if os.path.isfile(flagpath) and not unflag:
            flag = np.fromfile(flagpath, dtype='bool').reshape((nt, nf))
        else:
            flag = np.zeros((nt, nf), dtype='bool')

        if time is not None:
            flag = flag[time]
        elif frequency is not None:
            flag = flag[:,frequency]

        for i,u in enumerate(ubls):
            p = p + 1
            plt.subplot(len(ubls), len(datafiles), p)

            if time is not None:
                plotdata = data[time, i]
            elif frequency is not None:
                plotdata = data[:, i, frequency]
            else:
                plotdata = data[:, i]

            if not force_range:
                ranges = [np.percentile(plotdata, 5), np.percentile(plotdata, 95)]

            plotdata[flag] = np.nan
            if len(plotdata.shape) == 2:
                plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
                plt.colorbar()
            else:
                plt.plot(plotdata)
                plt.ylim(ranges)
            plt.title('%s, baseline [%.1f, %.1f, %.1f], %s'%(os.path.basename(datafile), info['ubl'][u][0], info['ubl'][u][1], info['ubl'][u][2], mode))
    if oppath is not None:
        plt.savefig(oppath, bbox_inches='tight')
    if not suppress:
        plt.show()
    else:
        plt.close()


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

        data.shape = (nt, (nf + nprefix))
        data = data[:, nprefix:]

        flagpath = datafile.replace(plottype, 'omniflag')
        if os.path.isfile(flagpath) and not unflag:
            flag = np.fromfile(flagpath, dtype='bool').reshape((nt, nf))
        else:
            flag = np.zeros((nt, nf), dtype='bool')

        if time is not None:
            plotdata = np.copy(data[time])
            flag = np.copy(flag[time])
        elif frequency is not None:
            plotdata = np.copy(data[:, frequency])
            flag = np.copy(flag[:, frequency])
        else:
            plotdata = np.copy(data)
        p = p + 1
        plt.subplot(1, len(datafiles), p)
        if not force_range:
            ranges = [np.percentile(plotdata, 5), np.percentile(plotdata, 95)]
        else:
            ranges[0] = max(ranges[0], 0)
        plotdata[flag] = np.nan
        if len(plotdata.shape) == 2:
            plt.imshow(plotdata, vmin = ranges[0], vmax = ranges[1], interpolation = 'nearest')
            plt.colorbar()
        else:
            plt.plot(plotdata)
            plt.ylim(ranges)
        plt.title('%s, chi^2'%os.path.basename(datafile))
    if oppath is not None:
        plt.savefig(oppath, bbox_inches='tight')
    if not suppress:
        plt.show()
    else:
        plt.close()

######################################################################
############## omniview uv ###################################
######################################################################

if plottype == 'uv':
    plotdata = omni.pick_slice_uvs(datafiles, pol, time, frequency)
    #try:
        #totalVisibilityId = info['totalVisibilityId']
    #except KeyError:
        #nant = int(np.ceil((len(plotdata)*2)**.5))
        #totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(nant)])
    omni.omniview(plotdata, info, oppath = oppath, suppress= suppress)
