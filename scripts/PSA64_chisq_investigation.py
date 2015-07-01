import numpy as np
import matplotlib.pyplot as plt
import struct, ephem, glob
import aipy as ap
from scipy import interpolate
import scipy.ndimage.filters as sfil
import omnical.calibration_omni as omni
julDelta = 2415020.

uv = ap.miriad.UV('/data4/paper/2012EoR/psa_live/psa6375/zen.2456375.59142.uvcRREcAC')
sa = ephem.Observer()
sa.lon = uv['longitu']
sa.lat = uv['latitud']
sun = ephem.Sun()
del(uv)

tag = 'feb2015'#'march2015v6'#'feb2015'#
script_ver = '2'
#get data for all freq and lst, only care abt first time slice of each file
NT = 100
NF = 203
CEILING = 2.e8
MIN_FILTER_SIZE = (1, 1)
HEAVY_FILTER_FREQ_RANGE = [60, 150]
HEAVY_MIN_FILTER_SIZE = (1,1)#(1, 15)
MEDIAN_FILTER_SIZE = (6,3)
for POL in ['xx', 'yy']:
    chisqfiles=sorted(glob.glob('/data4/paper/hz2ug/2015PSA64/zen*%s*%s.omnichisq'%(tag, POL)))

    data = np.zeros((NT, NF)) + CEILING

    print "WARNING! assuming same dof for all chisqs"
    dof = omni.read_redundantinfo(chisqfiles[0].replace('.omnichisq','.binfo'), DoF_only=True)
    for i, chisqfile in enumerate(chisqfiles):
        #print '.',
        rawd = np.fromfile(chisqfile,dtype='float32')
        
        sa.date = struct.unpack('d', struct.pack('ff', *rawd[:2].tolist()))[0] - julDelta
        lst = sa.sidereal_time()
        sun.compute(sa)
        if sun.alt < -0.1:
            nf = rawd[2]
            nt = len(rawd)/(3 + nf)
            rawd.shape = (nt, 3 + nf)
            
            flag = np.fromfile(chisqfile.replace('.omnichisq','.omniflag'),dtype='bool').reshape((nt,nf))[:, ::int(np.floor(nf/NF))][:,:NF].all(axis=0)

            d = np.min(rawd[:, 3::int(np.floor(nf/NF))][:, :NF], axis = 0)/dof
            #print flag.dtype, (d==0).dtype
            flag = flag|(d==0)
            update_t = int(np.floor(lst/(2*np.pi)*NT))
            data[update_t, ~flag] = np.min([data[update_t, ~flag], d[~flag]], axis=0)


    #model = np.outer(np.median(data, axis = 1), np.median(data, axis = 0))
    no_model = (data==CEILING)
    no_model_freq = (np.sum(no_model, axis = 0) > NT/3)
    no_model = no_model | no_model_freq[None, :]
    data[no_model] = CEILING
    model = sfil.minimum_filter(data, size = MIN_FILTER_SIZE)
    model[:, HEAVY_FILTER_FREQ_RANGE[0]:HEAVY_FILTER_FREQ_RANGE[1]] = sfil.minimum_filter(data, size = HEAVY_MIN_FILTER_SIZE)[:, HEAVY_FILTER_FREQ_RANGE[0]:HEAVY_FILTER_FREQ_RANGE[1]]
    model = sfil.median_filter(model, size = MEDIAN_FILTER_SIZE)
    
    model = model * np.median(data[~no_model]/model[~no_model])
    plt.subplot('211');plt.imshow(np.log10(data), vmin = 4, vmax = 5);
    plt.subplot('212');plt.imshow(np.log10(model), vmin = 4, vmax = 5);
    opmodel = np.zeros((len(model), NF + 3), dtype='float32')

    for i in range(NT):
        opmodel[i, 0:2] = struct.unpack('ff', struct.pack('d', i * 2*np.pi/NT))
    opmodel[:, 2] = NF
    opmodel[:, 3:] = model
    opmodel.tofile('/data4/paper/hz2ug/2015PSA64/PSA64_chisq_model_%s_v%s_%s.omnichisq'%(tag, script_ver, POL))
    plt.show()

    #plotting for one frequency
    data = {}
    for j,wantbin in enumerate([40, 82, 110, 165]):
        data[wantbin] = np.zeros(len(chisqfiles), dtype='float32')
        lst = np.zeros(len(chisqfiles), dtype='float32')
        sundown = np.zeros(len(chisqfiles), dtype='bool')
        for i, chisqfile in enumerate(chisqfiles):
            rawd = np.fromfile(chisqfile,dtype='float32')
            nf = rawd[2]
            flag = np.fromfile(chisqfile.replace('.omnichisq','.omniflag'),dtype='bool')
            
            sa.date = struct.unpack('d', struct.pack('ff', *rawd[:2].tolist()))[0] - julDelta
            lst[i] = sa.sidereal_time()
            sun.compute(sa)
            sundown[i] = sun.alt < -0.1
            d = rawd[wantbin + 3::nf+3]/dof
            try:
                data[wantbin][i] = np.min(d[~flag[wantbin::nf]])
                if sundown[i] and 0 < lst[i] and lst[i] < 2 and abs(data[wantbin][i] - 8e4) < 2e4:
                    print chisqfile
            except:
                data[wantbin][i] = np.nan

        plt.subplot('14%i'%(j+1));plt.scatter(lst[sundown], data[wantbin][sundown]);plt.scatter(np.arange(0,2*np.pi,2*np.pi/NT), model[:,wantbin], color='green');plt.scatter(np.arange(0,2*np.pi,2*np.pi/NT), 1.5 * model[:,wantbin], color='cyan');plt.ylim([0,1e5]);plt.title(wantbin)

    plt.show()
exit()

###Zaki's noise model##
nchan = 203
tint = 43
f = np.load('/data2/home/hz2ug/PSA64_plot_for_zaki/tsys_model_jy_corrected.npz')
noise = np.array(f['tsys_jy'] **2) / (tint*100./nchan*1e6) * 1077

noise_model = interpolate.interp2d(np.arange(0, nchan), f['lsts']/12.*np.pi, noise)

plt.subplot('121')
plt.scatter(lst[sundown], data[sundown]);plt.ylim([0,1e8])

plt.subplot('122')
plt.scatter(f['lsts']/12.*np.pi, noise[:,wantbin]);plt.ylim([0,1e7])

plt.show()

