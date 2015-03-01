import numpy as np
import matplotlib.pyplot as plt
import struct, ephem, glob
import aipy as ap
from scipy import interpolate
import scipy.ndimage.filters as sfil
julDelta = 2415020.

chisqfiles=sorted(glob.glob('/data4/paper/hz2ug/2015PSA64/*xx.omnichisq'))



uv = ap.miriad.UV('/data4/paper/2012EoR/psa_live/psa6375/zen.2456375.59142.uvcRREcAC')
sa = ephem.Observer()
sa.lon = uv['longitu']
sa.lat = uv['latitud']
sun = ephem.Sun()
del(uv)

#get data for all freq and lst, only care abt first time slice of each file
NT = 100
NF = 203
data = np.zeros((NT, NF)) + 1.e9
for i, chisqfile in enumerate(chisqfiles):
    rawd = np.fromfile(chisqfile,dtype='float32')
    
    sa.date = struct.unpack('d', struct.pack('ff', *rawd[:2].tolist()))[0] - julDelta
    lst = sa.sidereal_time()
    sun.compute(sa)
    if sun.alt < -0.1:
        nf = rawd[2]
        d = rawd[3:3+nf:int(np.floor(nf/NF))][:NF]
        flag = np.fromfile(chisqfile.replace('.omnichisq','.omniflag'),dtype='bool')[:nf:int(np.floor(nf/NF))][:NF]
        #print flag.dtype, (d==0).dtype
        flag = flag|(d==0)
        update_t = int(np.floor(lst/(2*np.pi)*NT))
        data[update_t, ~flag] = np.min([data[update_t, ~flag], d[~flag]], axis=0)


#model = np.outer(np.median(data, axis = 1), np.median(data, axis = 0))
model = sfil.minimum_filter(data, size = 5)

model = model * np.median(data/model)
plt.subplot('211');plt.imshow(np.log10(data), vmin = 7, vmax = 9);
plt.subplot('212');plt.imshow(np.log10(model), vmin = 7, vmax = 9);
model.tofile('/data4/paper/hz2ug/2015PSA64/chisq_model_t%i_f%i.bin'%(NT,NF))
plt.show()

#plotting for one frequency
for j,wantbin in enumerate([40, 80, 110, 165]):
    data = np.zeros(len(chisqfiles), dtype='float32')
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
        d = rawd[wantbin + 3::nf+3]
        try:
            data[i] = np.min(d[~flag[wantbin::nf]])
        except:
            data[i] = np.nan

    plt.subplot('14%i'%(j+1));plt.scatter(lst[sundown], data[sundown]);plt.ylim([0,1e8]);plt.scatter(np.arange(0,2*np.pi,2*np.pi/NT), 1.7 * model[:,wantbin]);plt.ylim([0,1e8]);plt.title(wantbin)

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

