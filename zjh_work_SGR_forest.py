
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os
import pandas as pd
from daily_search.tool import overlap_bins
from daily_search.perception import Event
from daily_search.background import get_background_f,TD_baseline
import h5py

marker = 'tg_20200427_18330585'

datatop = '/home/laojin/my_work/' + marker + '/'
savedir = '/home/laojin/my_work/' + marker + '/'

if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)


detectors = ['b1','n9','n6','n0']
binsize = 0.1
stepsize = 0.05
bins = overlap_bins([-15,30],binsize,stepsize)
lc_binsize = 0.064
bins_lc = np.arange(-17,31,lc_binsize)

for dete in detectors:
	link = datatop + marker +'_' + dete + '.fits'
	hl = fits.open(link)
	trigtime = hl[0].header['TRIGTIME']

	time = hl[2].data.field(0)
	ch = hl[2].data.field(1)
	ch_n = hl[1].data.field(0)
	e_ch = np.vstack([ch_n,ch_n]).T
	e1 = hl[1].data.field(1)
	e2 = hl[1].data.field(2)
	t = time - trigtime
	ch_E = pd.DataFrame({'CHANNEL':ch_n,'E_MIN':e1,'E_MAX':e2})
	lightcurve = Event(ch_E,t,ch)

	overlap_lc_t,overlap_lc_rates = lightcurve.ovelap_lightcurve(bins,energy_band=e_ch,channel= True)
	lc_t,lc_rates = lightcurve(bins = bins_lc,energy_band=e_ch,channel=True)
	bs_list = []
	dete_save = savedir + 'A_lc_ch_'+dete + '/'
	if os.path.exists(dete_save)  ==False:
		os.makedirs(dete_save)
	for index,lc_ratesi in enumerate(lc_rates):
		bs_f = get_background_f(lc_t,lc_ratesi)
		bs = bs_f(overlap_lc_t)
		bs_list.append(bs)

		plt.title('lc of channel '+str(index))
		plt.plot(overlap_lc_t,overlap_lc_rates[index],color = 'k',label = 'overlap lc')
		plt.plot(overlap_lc_t,bs,color = 'orange',label = 'overlap bs')
		plt.axvline(x = 0,color = 'r',label = 'trig time')
		plt.legend()
		plt.xlabel('time (s)')
		plt.ylabel('connet rate')
		plt.savefig(dete_save + 'channel_'+ str(index) + '.png')
		plt.close()

	bs_arr = np.vstack(bs_list)
	spectrum_save = savedir + 'Z_pha_'+ dete +'.hdf5'
	if os.path.exists(spectrum_save):
		os.remove(spectrum_save)
	spc = overlap_lc_rates.T

	all_rate = spc.sum(axis=1)
	bs = TD_baseline(overlap_lc_t,all_rate)
	SNR = (all_rate-bs)/np.sqrt(bs/binsize)

	spc_bs = bs_arr.T
	f1 = h5py.File(spectrum_save,'w')
	time_ = f1.create_dataset('time',overlap_lc_t.shape,dtype = np.float)
	time_[...] = overlap_lc_t
	SNR_ = f1.create_dataset('SNR',SNR.shape,dtype = np.float)
	SNR_[...] = SNR
	spc_ = f1.create_dataset('spc',spc.shape,dtype = np.float)
	spc_[...] = spc
	spc_bs_ = f1.create_dataset('spc_bs',spc_bs.shape,dtype = np.float)
	spc_bs_[...] = spc_bs
	f1.close()


