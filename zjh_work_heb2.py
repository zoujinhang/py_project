import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from Calculate_duration import get_txx
from Calculate_duration import Plot

savedir = '/home/laojin/my_work/heb/result/'
hebdatalink = '/home/laojin/my_work/heb/data/HEB190530429_HE-Evt_detc12.fits'
ttedatalink = '/home/laojin/trigdata/bn190530430/glg_tte_n1_bn190530430_v00.fit'

binsize = 1
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl = fits.open(hebdatalink)
heb_time = hl[1].data.field(0)
heb_trigtime = hl[0].header['TRIGTT']
heb_t = heb_time - heb_trigtime

heb_index = np.where((heb_t>=-60)&(heb_t<=110))[0]
heb_t = heb_t[heb_index]

tte = fits.open(ttedatalink)
trigtime = tte[0].header['TRIGTIME']
time = tte[2].data.field(0)
t = time - trigtime
tte_index = np.where((t>=-60)&(t<=110))[0]
t = t[tte_index]

heb_edges = np.arange(-50,100+binsize,binsize)
heb_n,heb_edges = np.histogram(heb_t,bins = heb_edges)
tte_n,tte_edges = np.histogram(t,bins=heb_edges)
heb_c = (heb_edges[1:]+heb_edges[:-1])*0.5
heb_rate = heb_n/binsize
tte_rate = tte_n/binsize
plt.figure()
plt.subplot(2,1,1)
plt.plot(heb_c,heb_rate)
plt.subplot(2,1,2)
plt.plot(heb_c,tte_rate)
plt.savefig(savedir+ 'A_HEB190530429_lightcurve_.png')
plt.close()

heb_results = get_txx(heb_t,binsize = binsize,time_unified=False,hardness = 100,it = 1000,block_n = 40,prior=12)
plt.title('HEB190530429')
myplt = Plot(heb_results)
myplt.plot_light_curve()
plt.xlim(-50,70)
plt.savefig(savedir + 'A_HEB190530429_lightcurve.png')
plt.close()

if heb_results['good']:
	for i in range(len(heb_results['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_HEB190530429_distribution_'+str(i)+'.png')
		plt.close()
plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
plt.title('HEB190530429')
myplt.plot_Txx1('90')
plt.xlim(-50,70)
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim(-50,70)
plt.savefig(savedir + 'B_HEB190530429_txx.png')
plt.close()

tte_results = get_txx(t,binsize = binsize,time_unified=False,hardness = 100,it = 1000,block_n = 40,prior=12)
plt.title('bn190530430')
myplt = Plot(tte_results)
myplt.plot_light_curve()
plt.xlim(-50,70)
plt.savefig(savedir + 'A_bn190530430_lightcurve.png')
plt.close()

if tte_results['good']:
	for i in range(len(tte_results['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_bn190530430_distribution_'+str(i)+'.png')
		plt.close()
plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
plt.title('bn190530430')
myplt.plot_Txx1('90')
plt.xlim(-50,70)
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim(-50,70)
plt.savefig(savedir + 'B_bn190530430_txx.png')
plt.close()

myplt.plot_normallization()
plt.savefig(savedir+'C_bn190530430_normall.png')
plt.close()

myplt.plot_ACC()
plt.savefig(savedir + 'D_bn190530430_ACC.png')
plt.close()


