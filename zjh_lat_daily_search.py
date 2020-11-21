
from Data_analysis import TD_baseline,TD_bs
import Data_analysis.file as myfile
import numpy as np
import pandas as pd
import os
from Fermi_tool.daily import Database
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
from Data_analysis.geometry import Geometry
from Data_analysis import ch_to_energy
from Fermi_tool.lag import get_band



topdir = '/home/laojin/daily_data/'
#topdir = '/media/laojin/Elements/daily/'
#topdir = '/media/laojin/TOSHIBA_EXT/daily/'
#timestart = '2020-04-28T00:00:00'
#timestop = '2020-04-28T00:59:00'

timestart = '2020-10-01T05:57:00'
timestop = '2020-10-01T05:59:00'

savedir = '/home/laojin/my_lat/daily_search/'
binsize = 0.064

e_band = get_band([8,800],5,ovelap=0)
print('e_band',e_band)


if os.path.exists(savedir) == False:
	os.makedirs(savedir)


detector_list = ['n0' ,'n1' ,'n2' ,'n3' ,'n4' ,'n5' ,'n6' ,'n7' ,'n8' ,'n9' ,'na' ,'nb']

geometry = Geometry()
fermi = Database(topdir)
data = fermi.get_detector_data(timestart ,timestop)

t_all = np.array([])

for ni_index in detector_list:

	ni = data[ni_index]
	events = ni['events']
	t = events['TIME'].values
	ch = events['PHA'].values
	bins_i = np.arange(t.min(), t.max() - 1, binsize)
	bin_c_i = 0.5 * (bins_i[1:] + bins_i[:-1])
	ch_n = ni['ch_E']['CHANNEL'].values
	e1 = ni['ch_E']['E_MIN'].values
	e2 = ni['ch_E']['E_MAX'].values

	t,e = ch_to_energy(t,ch,ch_n,e1,e2)

	plt.figure(figsize = (10,10))

	for index,band in enumerate(e_band):
		plt.subplot(len(e_band)+1, 1,len(e_band)+1-index)
		index_e = np.where((e>=band[0])&(e<band[-1]))[0]
		t_e = t[index_e]
		e_e = e[index_e]

		bn_n_i_e = np.histogram(t_e,bins = bins_i)[0]
		rate_i_e = bn_n_i_e/binsize
		cs2_i_e,bs2_i_e,sigma_i_e = TD_bs(bin_c_i,rate_i_e,lambda_=5000,sigma=True)
		plt.plot(bin_c_i-t.min(),rate_i_e,label = str(band))
		plt.plot(bin_c_i-t.min(),bs2_i_e,label = 'TD_bs')
		plt.plot(bin_c_i-t.min(),bs2_i_e+3*sigma_i_e,label = 'sigma')
		plt.legend()

	bn_n_i = np.histogram(t,bins = bins_i)[0]
	rate_i = bn_n_i/binsize

	bin_c_i,cs_i,bs_i = TD_baseline(bin_c_i,rate_i)
	cs2_i,bs2_i,sigma_i = TD_bs(bin_c_i,rate_i,lambda_=5000,sigma=True)

	plt.subplot(len(e_band)+1,1,1)
	plt.plot(bin_c_i-t.min(),rate_i)
	plt.plot(bin_c_i-t.min(),bs_i)
	plt.plot(bin_c_i-t.min(),bs2_i,label = 'TD_bs')
	plt.plot(bin_c_i-t.min(),bs2_i+3*sigma_i,label = 'sigma')
	plt.legend()
	plt.savefig(savedir + 'B_lightcurve_'+ni_index+'.png')
	plt.close()
	t_all = np.concatenate((t_all,t))

t_all = np.sort(t_all)

bins = np.arange(t_all.min(),t_all.max()-1,binsize)

bn_n = np.histogram(t_all,bins = bins)[0]
rate = bn_n/binsize
bin_c = 0.5*(bins[1:]+bins[:-1])

bin_c,cs,bs = TD_baseline(bin_c,rate)
cs2,bs2,sigma_ = TD_bs(bin_c,rate,lambda_=5000,sigma=True)
plt.plot(bin_c-t_all.min(),rate,label = 'lightcurve')
plt.plot(bin_c-t_all.min(),bs,label = 'TD_baseline')
plt.plot(bin_c-t_all.min(),bs2,label = 'TD_bs')
plt.plot(bin_c-t_all.min(),bs2+3*sigma_,label = 'sigma')
#plt.xlim(1060,1065)
plt.legend()
plt.savefig(savedir + 'A_lightcurve111.png')
plt.close()













