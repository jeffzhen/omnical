#!/usr/bin/env python

import aipy as ap
import numpy as np
import commands, os, time, math, ephem, optparse, sys, warnings
import omnical.calibration_omni as omni
import cPickle as pickle
import scipy.signal as ss
import scipy.ndimage.filters as sfil
from scipy import interpolate
FILENAME = "omnical_PSA128.py"
print "#Omnical Version %s#"%omni.__version__
PI = np.pi
TPI = 2 * np.pi
TIMER = omni.Timer()
######################################################################
##############Config parameters###################################
######################################################################
o = optparse.OptionParser()

#ap.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-t', '--tag', action = 'store', default = 'PSA128', help = 'tag name of this calibration')
o.add_option('-C', '--cal', action = 'store', default = None, help = 'calfile for processing uv files')
o.add_option('-d', '--datatag', action = 'store', default = None, help = 'tag name of this data set')
o.add_option('-o', '--outputpath', action = 'store', default = ".", help = 'output folder')
o.add_option('-f', '--overwrite', action = 'store_true', help = 'whether to overwrite if the new uv files already exists')
o.add_option('--skip_sun', action = 'store_true', help = 'whether to calibrate data set with sun up.')
o.add_option('--mem', action = 'store', type = 'float', default = 4e9, help = 'Amount of initial memory to reserve when parsing uv files in number of bytes.')
o.add_option('--compress_bin', action = 'store', type = 'int', default = 19, help = 'Compress each uvfile (not including first and last) to the specified number of time bins.')

opts,args = o.parse_args(sys.argv[1:])
calibrate_sun = not opts.skip_sun
overwrite_uvs = opts.overwrite
ano = opts.tag##This is the file name difference for final calibration parameter result file. Result will be saved in miriadextract_xx_ano.omnical
dataano = opts.datatag#ano for existing data and lst.dat
oppath = os.path.expanduser(opts.outputpath)
delay_compression = opts.compress_bin
uvfiles = [os.path.expanduser(arg) for arg in args]

for uvf in uvfiles:
    if not os.path.isdir(uvf):
        uvfiles.remove(uvf)
        print "WARNING: uv file path %s does not exist!"%uvf
if len(uvfiles) < 3:
    raise Exception("ERROR: Not enough valid uv files detected (%i < 3) in input. Exiting!"%len(uvfiles))


init_mem = opts.mem


if dataano is None:
    dataano = ''
    for i, uvf in enumerate(uvfiles[1:-1]):
        if i!= 0:
            dataano = dataano + '_'
        while os.path.basename(uvf) == '' and len(uvf) > 0:
            uvf = uvf[:-1]
        dataano = dataano + os.path.basename(uvf)

print "Reading calfile %s"%opts.cal,
sys.stdout.flush()
aa = ap.cal.get_aa(opts.cal, np.array([.15]))
print "Done"
sys.stdout.flush()

######################################################################
######################################################################
######################################################################

########Massage user parameters###################################
oppath += '/'

####get some info from the first uvfile   ################
print "Getting some basic info from %s"%uvfiles[0],
sys.stdout.flush()

sa = ephem.Observer()
sa.pressure = 0
uv=ap.miriad.UV(uvfiles[0])
nfreq = uv.nchan;
nant = uv['nants']
sa.lon = aa.lon
sa.lat = aa.lat
startfreq = uv['sfreq']#GHz
dfreq = uv['sdf']#GHz
wantpols = {}
for pol in omni.get_uv_pols(uv):
    wantpols[pol] = ap.miriad.str2pol[pol]
del(uv)
print "Done."
sys.stdout.flush()


#######import data file###########################
print FILENAME + " MSG:",  len(uvfiles), "uv files to be processed for " + ano
sys.stdout.flush()
data, jd, timing, lst, rawflag = omni.importuvs(uvfiles, wantpols, totalVisibilityId = np.concatenate([[[i,j] for i in range(j + 1)] for j in range(len(aa))]), timingTolerance=100, init_mem=init_mem, lat = sa.lat, lon=sa.lon)#, nTotalAntenna = len(aa))
if delay_compression > len(jd) / len(uvfiles):
    raise ValueError("Desired time bins %i is larger than input stamp count per uvfile."%len(jd))

data[rawflag] = 0
print FILENAME + " MSG:",  len(jd), "slices read. data shape: ", data.shape
sys.stdout.flush()


########compress###########################
compr_shape = list(data.shape)
compr_shape[1] = len(uvfiles) * delay_compression
compressed_data = np.zeros(compr_shape, dtype=data.dtype)
compr_flags = np.zeros(compr_shape, dtype='bool')
compr_error = np.zeros((compr_shape[0], compr_shape[2], compr_shape[-1]), dtype='float32')
for p, pol in enumerate(wantpols.keys()):
    print FILENAME + " MSG: starting compression on %s %s"%(dataano, pol),
    sys.stdout.flush()
    timer = time.time()

    for f in range(data.shape[2]):
        print '.',
        sys.stdout.flush()
        flag = rawflag[p, :, f].all(axis = -1)
        badbl = np.sum(rawflag[p, :, f], axis = 0) > np.sum(flag)
        compr_result = omni.deconvolve_spectra2(data[p, :, f], ~flag, (compr_shape[1]+1)/2, correction_weight=1e-6)
        compressed_data[p, :, f] = float(compr_shape[1])/nfreq * compr_result[0]
        compr_error[p,f] = compr_result[1]
        compr_error_bar = float(compr_shape[1])/nfreq * np.array([compr_result[3][i,i] for i in range(compr_shape[1])])**.5
        compr_flags[p, compr_error_bar > min(compr_error_bar) * 1.05, f] = True
        compr_flags[p, :, f, badbl] = True

    print "Done. %fmin"%(float(time.time()-timer)/60.)
    sys.stdout.flush()
#######save uvs
print FILENAME + " MSG: saving compressed uv files",
sys.stdout.flush()
if len(wantpols.keys()) == 1:
    compr_uvfile = oppath + '/' + dataano + "_%s.uvOEE"%pol
else:
    compr_uvfile = oppath + '/' + dataano + ".uvOEE"


uv=ap.miriad.UV(uvfiles[0])
jd=np.array(jd)
deljd = np.mean(jd[1:] - jd[:-1])
new_jd = np.arange(compr_shape[1])**2 * deljd / data.shape[1]
omni.exportuv(compr_uvfile, compressed_data[:,delay_compression:-delay_compression], compr_flags[:,delay_compression:-delay_compression], wantpols.values(), new_jd[delay_compression:-delay_compression], uv['inttime'], uv['sfreq'], uv['sdf'], sa.lat, sa.lon, overwrite=False, comment=ano)
del(uv)
print "Done."
sys.stdout.flush()
