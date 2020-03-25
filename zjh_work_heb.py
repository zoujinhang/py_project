import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from Calculate_duration import get_txx
from Calculate_duration import Plot

savedir = '/home/laojin/my_work/heb/result/'
hebdatalink = '/home/laojin/my_work/heb/data/HEB171223818_HE-Evt_detec09.fits'
ttedatalink = '/home/laojin/trigdata/bn171223818/glg_tte_n6_bn171223818_v00.fit'

binsize = 0.04


if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl =  fits.open(hebdatalink)
heb_time = hl[0].data.T[0]
heb_time  = np.sort(heb_time)
heb_time = heb_time-188681898.00

heb_index = np.where((heb_time>=-5)&(heb_time<=5.1))[0]
heb_time = heb_time[heb_index]

tte = fits.open(ttedatalink)
trigtime = tte[0].header['TRIGTIME']
time = tte[2].data.field(0)
t = time - trigtime
tte_index = np.where((t>=-10)&(t<=10))[0]
t = t[tte_index]

print(heb_time)
print(heb_time[-1]-heb_time[0])

heb_edges = np.arange(-2,2+binsize,binsize)
heb_n,heb_edges = np.histogram(heb_time,bins = heb_edges)

tte_n,tte_edges = np.histogram(t,bins=heb_edges)

heb_c = (heb_edges[1:]+heb_edges[:-1])*0.5
heb_rate = heb_n/binsize
tte_rate = tte_n/binsize

plt.figure()
plt.subplot(2,1,1)
plt.plot(heb_c,heb_rate)
plt.xlim(-1,1)
plt.subplot(2,1,2)
plt.plot(heb_c,tte_rate)
plt.xlim(-1,1)
plt.savefig(savedir+ 'A_heb_lightcurve.png')
plt.close()

heb_results = get_txx(heb_time,binsize = binsize,time_unified=False,hardness = 100,it = 1000,block_n = 40,prior=12)
plt.title('HEB171223818')
myplt = Plot(heb_results)
myplt.plot_light_curve()
plt.xlim(-1,1)
plt.savefig(savedir + 'A_lightcurve1.png')
plt.close()

if heb_results['good']:
	for i in range(len(heb_results['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_distribution1_'+str(i)+'.png')
		plt.close()
plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
plt.title('HEB171223818')
myplt.plot_Txx1('90')
plt.xlim([-1,1])
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim([-1,1])
plt.savefig(savedir + 'B_txx1.png')
plt.close()

myplt.plot_normallization()
plt.savefig(savedir+'C_normall1.png')
plt.close()

myplt.plot_ACC()
plt.savefig(savedir + 'D_ACC1.png')
plt.close()


tte_results = get_txx(t,binsize = binsize,time_unified=False,hardness = 100,it = 1000,block_n = 40,prior=12)
plt.title('bn171223818')
myplt = Plot(tte_results)
myplt.plot_light_curve()
plt.xlim(-1,1)
plt.savefig(savedir + 'A_lightcurve2.png')
plt.close()

if tte_results['good']:
	for i in range(len(tte_results['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_distribution2_'+str(i)+'.png')
		plt.close()
plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
plt.title('bn171223818')
myplt.plot_Txx1('90')
plt.xlim([-1,1])
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim([-1,1])
plt.savefig(savedir + 'B_txx2.png')
plt.close()

myplt.plot_normallization()
plt.savefig(savedir+'C_normall2.png')
plt.close()

myplt.plot_ACC()
plt.savefig(savedir + 'D_ACC2.png')
plt.close()








