import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from Calculate_duration import Plot,get_txx,save_result
from Data_analysis import ch_to_energy

datalink = '/media/laojin/Elements/trigdata/2020/bn200415367/glg_tte_n1_bn200415367_v00.fit'
savedir = '/home/laojin/my_lat/bn200415367_duration/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl = fits.open(datalink)

trigtime = hl[0].header['TRIGTIME']

time = hl[2].data.field(0)
ch = hl[2].data.field(1)
ch_n = hl[1].data.field(0)
e1 = hl[1].data.field(1)
e2 = hl[1].data.field(2)
t = time - trigtime
t, e = ch_to_energy(t, ch, ch_n, e1, e2)
e_lin = [[10,1000],[50,300],[10,350]]

for i in e_lin:
	e_index = np.where((e >= i[0]) & (e <= i[1]))[0]
	t = t[e_index]
	e = e[e_index]
	
	t_start = -0.5
	t_stop = 0.5
	result = get_txx(t,binsize = 0.064,time_edges=[t_start,t_stop],background_degree=7,sigma = 4,txx = 0.9,it = 300,p0 = 0.05,plot_check=savedir + 'Z_bn200415367_'+str(i[0])+'_'+str(i[1])+'_check.png',hardnss=100)
	save_result(result,savedir + 'C_bn200415367_'+str(i[0])+'_'+str(i[1])+'_T90.csv',float_format='%.3f')
	myplt = Plot(result)
	plt.title('GRB200415367 '+str(i[0])+'--'+str(i[1])+' kev')
	myplt.plot_light_curve(sigma=4.5)
	plt.xlim(-0.5,0.5)
	plt.savefig(savedir + 'A_bn200415367_'+str(i[0])+'_'+str(i[1])+'_lightcurve.png')
	plt.close()
	if result['good']:
		for ij in range(len(result['txx'])):
			plt.title('GRB200415367')
			myplt.plot_distribution('90',num = ij)
			plt.savefig(savedir + 'D_bn200415367_distribution_'+str(i[0])+'_'+str(i[1])+'_'+str(ij)+'.png')
			plt.close()
	#--------------------------------------------------------------------------
	#
	plt.figure(figsize = (10,10))
	plt.subplot(2,1,1)
	plt.title('GRB200415367 '+str(i[0])+'--'+str(i[1])+' kev')
	myplt.plot_Txx1('90')#
	plt.xlim(-0.5,0.5)
	plt.subplot(2,1,2)
	myplt.plot_Txx2('90')
	plt.xlim(-0.5,0.5)
	plt.savefig(savedir + 'B_bn200415367_'+str(i[0])+'_'+str(i[1])+'_txx.png')
	plt.close()






