import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from Calculate_duration import get_txx
from Calculate_duration import Plot


datalink = '/media/laojin/Elements/trigdata/2020/bn200219317/glg_tte_n3_bn200219317_v00.fit'

savedir = '/home/laojin/result/bn200219317_zbb/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hdu = fits.open(datalink)
trigtime = hdu[0].header['TRIGTIME']

ch_n = hdu[1].data.field(0)
e1 = hdu[1].data.field(1)
e2 = hdu[1].data.field(2)

time = hdu[2].data.field(0)
ch = hdu[2].data.field(1)
t = time - trigtime

e_up1 = 90
e_dw1 = 11

e_up2 = 127
e_dw2 = 7
binsize = 0.064

index_1 = np.where((ch >= e_dw1) & (ch <= e_up1))[0]
t_1 = t[index_1]

index_2 = np.where((ch >= e_dw2) & (ch <= e_up2))[0]
t_2 = t[index_2]

results1 = get_txx(t_1,binsize = binsize,block_time = 20,time_unified=False,hardness = 500*1/binsize,sigma = 1,it = 1000)
myplt = Plot(results1)
myplt.plot_light_curve()
plt.savefig(savedir + 'A_lightcurve1.png')
plt.close()

if results1['good']:
	for i in range(len(results1['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_distribution1_'+str(i)+'.png')
		plt.close()

plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
myplt.plot_Txx1('90')
plt.xlim([-5,5])
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim([-5,5])
plt.savefig(savedir + 'B_txx_15-350.png')
plt.close()

myplt.plot_normallization()
plt.savefig(savedir+'C_normall1.png')
plt.close()

myplt.plot_ACC()
plt.savefig(savedir + 'D_ACC1.png')
plt.close()

results2 = get_txx(t_2,binsize = binsize,block_time = 30,time_unified=False,hardness = 500*1/binsize,sigma = 1,it = 1000)
myplt = Plot(results2)
myplt.plot_light_curve()
plt.savefig(savedir + 'A_lightcurve2.png')
plt.close()

if results2['good']:
	for i in range(len(results2['txx'])):
		myplt.plot_distribution('90',num = i)
		plt.savefig(savedir + 'D_distribution2_'+str(i)+'.png')
		plt.close()

plt.figure(figsize = (10,10))
plt.subplot(2,1,1)
myplt.plot_Txx1('90')
plt.xlim([-5,5])
plt.subplot(2,1,2)
myplt.plot_Txx2('90')
plt.xlim([-5,5])
plt.savefig(savedir + 'B_txx_10-1000.png')
plt.close()

myplt.plot_normallization()
plt.savefig(savedir+'C_normall2.png')
plt.close()

myplt.plot_ACC()
plt.savefig(savedir + 'D_ACC2.png')
plt.close()



