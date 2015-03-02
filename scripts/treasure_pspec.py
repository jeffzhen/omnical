import numpy as np
import commands, os, time, math, ephem
import omnical.calibration_omni as omni
import omnical._omnical as _O
import optparse, sys
import matplotlib.pyplot as plt
import numpy.fft as fft
import scipy.signal as ss
import scipy.linalg as sla
import numpy.linalg as la
Cspeed = .299792458 #m/ns

treasure = omni.Treasure("/home/omniscope/data/PAPER/2015PSA64.treasure")
pol = 'xx'

#for ubl in treasure.ubls[pol]:
	#f=treasure.seal_name((pol, ubl))
	#if treasure.get_coin((pol, ubl)) is None:
		#print ubl,f
#exit()

#data = []
#for ubl in treasure.ubls[pol]:
	#c = treasure.get_coin_now((pol, ubl))
	#if c is not None:
		#plt.plot(abs(fft.fft(c.weighted_mean[10, 20:170])))
#plt.show()
#exit()

#for i, ubl in enumerate(treasure.ubls[pol][np.argsort(la.norm(treasure.ubls[pol],axis=1))[[4,24,44,64,84]]]):
	#c = treasure.get_coin_now((pol, ubl))
	#if c is not None:
		#plt.subplot('18%i'%(i+1))
		#plt.imshow(abs(fft.fft(c.weighted_mean[:30, 20:170], axis=1)), vmax=1e4)
#plt.show()
#exit()

ubl = treasure.ubls[pol][np.argsort(la.norm(treasure.ubls[pol],axis=1))[4]]#4
coin = treasure.get_coin_now((pol, ubl))
freqs = .1 + .1 * np.arange(coin.count.shape[1]) / coin.count.shape[1] #GHz
fbin_min = 20
fbin_max = 170
nf = fbin_max - fbin_min
roll_t = 150
tbin_min = 00
tbin_max = 1500#len(coin.count)
vertical_cut_thresh = 350#number of time slices flags above which all times will be flagged
delay_width = 1 / (freqs[fbin_max] - freqs[fbin_min]) #in nanoseconds
foreground_delay_thresh = la.norm(ubl) / Cspeed + delay_width * 2
foreground_deconvolve_delay_thresh = la.norm(ubl) / Cspeed + delay_width * 2
delays = delay_width * np.roll(np.arange(-np.floor((nf-1)/2.), np.ceil((nf+1)/2.)), -int(np.floor((nf-1)/2)))

print "ubl is %s, %.2fns, delay bin width is %.2fns."%(ubl, la.norm(ubl) / Cspeed, delay_width)

#plt.subplot('151')
#plt.imshow(coin.count, interpolation ='none'); plt.colorbar()
#plt.subplot('152')
#plt.imshow(abs(coin.mean), interpolation ='none',vmax=1e3); plt.colorbar()
#plt.subplot('153')
#plt.imshow(abs(coin.weighted_mean), interpolation ='none',vmax=1e3); plt.colorbar()
#plt.subplot('154')
#plt.imshow(coin.variance_re**.5, interpolation ='none',vmax=5); plt.colorbar()
#plt.subplot('155')
#plt.imshow(coin.weighted_variance**.5, interpolation ='none',vmax=5); plt.colorbar()
#plt.show()
#exit()


data = np.roll(coin.weighted_mean, roll_t, axis = 0)[tbin_min:tbin_max, fbin_min:fbin_max]

var = np.roll(coin.weighted_variance, roll_t, axis = 0)[tbin_min:tbin_max, fbin_min:fbin_max]


flag = np.isnan(data*var)|np.isinf(data*var)|(var > 15)
flag = flag | ((np.sum(flag[1200:1400], axis=0) > 50)[None,:])

#plt.plot(np.sum(flag, axis=0))
#plt.plot(np.sum(flag[1200:1400], axis=0))

#plt.show()
#exit()

#flag = flag|(var > 15)
#data = np.ones_like(data) * 1.e3 * 5e4 * np.arange(200.+fbin_min, 200+fbin_max)[None,:]**-2
#data[flag] = np.nan
#data = data - np.nanmean(data, axis = 0)[None,:]
data[flag] = 0
bhfilter = ss.blackmanharris(nf)
fdata = fft.fft(data * bhfilter[None, :], axis=1)
fwindow = fft.fft((~flag) * bhfilter[None, :], axis=1)
full_deconv_fdata = np.zeros_like(fdata)
#qaz = np.zeros_like(fdata)
#model_qaz = np.zeros_like(fdata)
model_fdata = np.zeros_like(fdata)

#i = 600
#m = sla.toeplitz(fwindow[i], np.roll(fwindow[i][::-1], 1)).astype('complex128')[:, abs(delays) <= foreground_deconvolve_delay_thresh]
#el, ec = sla.eigh(m.transpose().conjugate().dot(m))
#plt.plot(np.log10(el))
#plt.show();exit()

for i in range(len(fdata)):
	m = sla.toeplitz(fwindow[i], np.roll(fwindow[i][::-1], 1)).astype('complex128')[:, abs(delays) <= foreground_deconvolve_delay_thresh]
	mmi = la.inv(m.transpose().conjugate().dot(m) + np.identity(m.shape[1])*1e-2)
	deconv_fdata = mmi.dot(m.transpose().conjugate().dot(fdata[i]))
	full_deconv_fdata[i, abs(delays) <= foreground_deconvolve_delay_thresh] = deconv_fdata
	full_deconv_fdata[i, abs(delays) > foreground_delay_thresh] = 0
	model_fdata[i] = m.dot(full_deconv_fdata[i, abs(delays) <= foreground_deconvolve_delay_thresh])#m.dot(deconv_fdata)
model_data = fft.ifft(full_deconv_fdata,  axis=1) * (fbin_max - fbin_min)
residual= data-model_data
residual[flag] = np.nan


###r = int(np.floor((nf-1)/2))
###r2 = len(full_deconv_fdata[0, abs(delays) <= foreground_deconvolve_delay_thresh]) / 2
###plt.subplot('151')
###plt.imshow(np.abs(np.roll(fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=1e4);plt.colorbar();plt.title('data')
###plt.subplot('152')
###plt.imshow(np.abs(np.roll(full_deconv_fdata[:, abs(delays) <= foreground_deconvolve_delay_thresh], r2, axis=1)), aspect = 1/5.*r2/r, interpolation='none', vmax=np.percentile(abs(full_deconv_fdata[:, abs(delays) <= foreground_delay_thresh]), 85));plt.colorbar();plt.title('deconvolved data')
###plt.subplot('153')
###plt.imshow(np.abs(np.roll(model_fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=1e4);plt.colorbar();plt.title('convolved (deconvolved data inside horizon)')
###plt.subplot('154')
###plt.imshow(np.abs(np.roll(model_fdata-fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=5e2);plt.colorbar();plt.title('residual')
###plt.subplot('155')
###plt.imshow(np.abs(np.roll(model_qaz, r, axis=1)), aspect = 1/5., interpolation='none', vmax=5e2);plt.colorbar();plt.title('residual')
###plt.show();exit()


r = int(np.floor((nf-1)/2))
r2 = len(full_deconv_fdata[0, abs(delays) <= foreground_deconvolve_delay_thresh]) / 2
plt.subplot('251')
plt.imshow(np.abs(data), aspect = 1/5., interpolation='none', vmax=2e3);plt.colorbar();plt.title('original data')
plt.subplot('256')
plt.imshow(np.abs(np.roll(fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=1e4);plt.colorbar();plt.title('data')
plt.subplot('257')
plt.imshow(np.abs(np.roll(fwindow, r, axis=1)), aspect = 1/5., interpolation='none', vmax=10);plt.colorbar();plt.title('windows function')
plt.subplot('258')
plt.imshow(np.abs(np.roll(full_deconv_fdata[:, abs(delays) <= foreground_deconvolve_delay_thresh], r2, axis=1)), aspect = 1/5.*r2/r, interpolation='none', vmax=np.percentile(abs(full_deconv_fdata[:, abs(delays) <= foreground_delay_thresh]), 85));plt.colorbar();plt.title('deconvolved data')
plt.subplot('259')
plt.imshow(np.abs(np.roll(model_fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=1e4);plt.colorbar();plt.title('convolved (deconvolved data inside horizon)')
plt.subplot(2,5,10)
plt.imshow(np.abs(np.roll(model_fdata-fdata, r, axis=1)), aspect = 1/5., interpolation='none', vmax=5e2);plt.colorbar();plt.title('residual')
plt.subplot('252')
plt.imshow(np.abs(model_data), aspect = 1/5., interpolation='none', vmax=2e3);plt.colorbar();plt.title('deconvolved data')
plt.subplot('253')
plt.imshow(np.log10(np.abs(residual)), aspect = 1/5., interpolation='none', vmin = 0, vmax=3);plt.colorbar();plt.title('residual amp')
plt.subplot('254')
plt.imshow(np.angle(residual), aspect = 1/5., interpolation='none', vmin = -np.pi, vmax=np.pi);plt.colorbar();plt.title('residual phase')
plt.subplot('255')
plt.imshow(np.log10(np.roll(coin.variance_re, roll_t, axis = 0)[tbin_min:tbin_max, fbin_min:fbin_max]**.5), aspect = 1/5., interpolation='none', vmin = 0, vmax=3);plt.colorbar()

#print full_deconv_fdata[300]
plt.show()
######i = 0
######plt.subplot('711')
######plt.plot(freqs[fbin_min:fbin_max], np.real(data[i]))
######plt.plot(freqs[fbin_min:fbin_max], np.imag(data[i]))
######plt.subplot('712')
######plt.plot(delays, abs(fwindow[i]))
######plt.subplot('713')
######plt.plot(delays, abs(fdata[i]))
######plt.subplot('714')
######full_deconv_fdata = np.zeros_like(fdata[i])
######full_deconv_fdata[abs(delays) <= foreground_delay_thresh] = deconv_fdata
######plt.plot(delays, abs(full_deconv_fdata))
######plt.subplot('715')
######plt.plot(delays, abs(m.dot(deconv_fdata)))
######plt.subplot('716')
######plt.plot(delays, abs(fdata[i] - m.dot(deconv_fdata)))
######plt.subplot('717')
######plt.plot(freqs[fbin_min:fbin_max], np.real(fft.ifft(full_deconv_fdata)))
######plt.plot(freqs[fbin_min:fbin_max], np.imag(fft.ifft(full_deconv_fdata)))
######plt.show()
######exit()




