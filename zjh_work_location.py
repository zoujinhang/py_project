import numpy as np
from daily_search import Database
from daily_search.file import readcol
from daily_search.background import TD_baseline
from daily_search.perception import Event
from daily_search.satellite import Geometry,Locate,Sky_map
import matplotlib
matplotlib.use('Agg')
import os,sys
import pandas as pd
import matplotlib.pyplot as plt

#para = sys.argv
#n_para = len(para)
#para_file = para[1]
#if n_para >2:
#	save_make = sys.argv[2]
#else:
#	save_make = None

savedir = '/home/laojin/my_lat/location_bb/'
#savetop = '/savetop/'
#datafile = '/para_file/' + para_file
#datatop = '/data-share/SSS_SHARE/DATA/SHAO_DOWN/data/'
datatop = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'

#if save_make is not None:
#	savedir = savetop + '/' + save_make + '/'
#else:
#	savedir = savetop

#utc_list,t1_list,t2_list = readcol(datafile)

utc_list = ['2020-04-27T21:00:00']
t1_list = [-5]
t2_list = [5]
for i_ in range(len(utc_list)):

	savedir_i = savedir + 'time_range_'+ str(i_) + '/'
	t1 = t1_list[i_]
	t2 = t2_list[i_]
	utc_trig = utc_list[i_]

	if os.path.exists(savedir_i) ==False:
		os.makedirs(savedir_i)

	t_eva = 0.5*(t2+t1)
	time_window_duration = (t2-t1)/0.5
	if time_window_duration<10.:
		time_window_duration = 10.

	t_w_1 = t_eva - 0.5*time_window_duration
	t_w_2 = t_eva + 0.5*time_window_duration


	fermi_data = Database(datatop)
	clock = fermi_data.clock

	t_trg_met = clock.utc_to_met(utc_trig)
	t_w_1_met = t_trg_met + t_w_1
	t_w_2_met = t_trg_met + t_w_2

	t1_met = t_trg_met + t1
	t2_met = t_trg_met + t2

	t_start = clock.met_to_utc(t_w_1_met)
	t_stop = clock.met_to_utc(t_w_2_met)

	data = fermi_data.get_detector_data(t_start,t_stop)
	pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)

	binsize_lightcurve = 0.05
	data_name = np.array(list(data.keys()))
	name = data_name[:12]
	geometry = Geometry(pd_pos_data)
	location = Locate(geometry)

	energy_band = [[8,105],[8,32],[33,84]]
	cenergies = [[5,50],[44,300]]

	bins_lightcurve = np.arange(t_w_1_met,t_w_2_met,binsize_lightcurve)
	bins_lightcurve_c = 0.5*(bins_lightcurve[1:]+bins_lightcurve[:-1])
	SNR_list = []
	lc_list = []
	bs_list = []

	for detei in name:
		ni = data[detei]
		ch_E = ni['ch_E']
		t_ = ni['events']['TIME'].values
		ch = ni['events']['PHA'].values
		lightcurve = Event(ch_E,t_,ch)
		detector_SNR_list = []
		detector_lc_list = []
		detector_bs_list = []
		for e_band in energy_band:
			t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=e_band,channel = True)
			bs = TD_baseline(t_cc,rate_c)
			scale = np.sqrt(bs/binsize_lightcurve)
			detector_SNR_list.append((rate_c-bs)/scale)
			detector_bs_list.append(bs)
			detector_lc_list.append(rate_c)

		SNR_list.append(detector_SNR_list)
		lc_list.append(detector_lc_list)
		bs_list.append(detector_bs_list)

	SNR_arr = []
	lc_arr = []
	bs_arr = []

	for i in range(len(energy_band)):
		e_band_SNR = []
		e_band_bs = []
		e_band_lc = []

		for j in range(len(name)):
			e_band_SNR.append(SNR_list[j][i])
			e_band_bs.append(bs_list[j][i])
			e_band_lc.append(lc_list[j][i])

		e_band_SNR = np.vstack(e_band_SNR).T
		e_band_bs = np.vstack(e_band_bs).T
		e_band_lc = np.vstack(e_band_lc).T

		SNR_arr.append(e_band_SNR)
		lc_arr.append(e_band_lc)
		bs_arr.append(e_band_bs)

	location5_50 = []
	location50_300 = []


	index_ = np.where((bins_lightcurve_c>=t1_met-1)&(bins_lightcurve_c<=t2_met+1))[0]

	index_list5_50 = []
	index_snr5_50 = None
	index_list50_300 = []
	index_snr50_300 = None

	for i in index_:
		index_snr = np.where(SNR_arr[1][i]>=5)[0]
		if len(index_snr)>=3 and len(index_list5_50)*binsize_lightcurve <= 2:
			index_list5_50.append(i)
			if index_snr5_50 is None:
				index_snr5_50 = index_snr
			else:
				if len(index_snr)>len(index_snr5_50):
					index_snr5_50=index_snr

		index_snr = np.where(SNR_arr[2][i]>=4.5)[0]
		if len(index_snr)>=3 and len(index_list50_300)*binsize_lightcurve <= 2:
			index_list50_300.append(i)
			if index_snr50_300 is None:
				index_snr50_300 = index_snr
			else:
				if len(index_snr)>len(index_snr50_300):
					index_snr50_300=index_snr

	if len(index_list5_50)>0:
		during = len(index_list5_50) * binsize_lightcurve
		if during>0:
			lc_5_50 = lc_arr[1][index_list5_50]
			m_rate_5_50 = lc_5_50.mean(axis = 0)
			bs_5_50 = bs_arr[1][index_list5_50]
			m_bs_5_50 = bs_5_50.mean(axis = 0)
			t = (bins_lightcurve_c[index_list5_50]).mean()
			detectorlist = name[index_snr5_50]
			location1 = location.locate(t,m_rate_5_50,m_bs_5_50,0,cenergies[0],detector_list=detectorlist,during=during)
			location5_50.append(location1)

		else:

			location5_50.append(None)

	else:
		location5_50.append(None)

	if len(index_list50_300)>0:

		during = len(index_list50_300) * binsize_lightcurve
		if during > 0:
			lc_50_300 = lc_arr[2][index_list50_300]
			m_rate_50_300 = lc_50_300.mean(axis = 0)
			bs_50_300 = bs_arr[2][index_list50_300]
			m_bs_50_300 = bs_50_300.mean(axis = 0)
			t = (bins_lightcurve_c[index_list50_300]).mean()
			detectorlist = name[index_snr50_300]
			location2 = location.locate(t,m_rate_50_300,m_bs_50_300,1,cenergies[1],detector_list=detectorlist,during=during)
			location50_300.append(location2)
		else:
			location50_300.append(None)
	else:
		location50_300.append(None)

	smp = Sky_map(figsize=(10,10))
	smp.add_subplot(2,1,1)
	smp.title('Energy band 5-50 kev')
	if location5_50[0] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location5_50[0]
		d_chi2 = chi2-chi2.min()
		smp.tricontour(ra_rcat,dec_rcat,d_chi2,[0,2.3,4.61,9.21])
		c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
		smp.scatter(ra_rcat, dec_rcat, marker=',', c=c, s=5)
		smp.plot(ra, dec, '*', markersize=10, color='orange')

	smp.plot_earth(t_trg_met,geometry)
	smp.plot_detector(t_trg_met,geometry)
	smp.plot_galactic_plane()
	smp.plot_sum(t_trg_met,geometry)
	smp.plot_moon(t_trg_met,geometry)
	smp.plot_continue_source()

	smp.add_subplot(2,1,2)
	smp.title('Energy band 50-300 kev')
	if location50_300[0] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location50_300[0]
		d_chi2 = chi2-chi2.min()
		smp.tricontour(ra_rcat,dec_rcat,d_chi2,[0,2.3,4.61,9.21])
		c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
		smp.scatter(ra_rcat, dec_rcat, marker=',', c=c, s=5)
		smp.plot(ra, dec, '*', markersize=10, color='orange')
	smp.plot_earth(t_trg_met,geometry)
	smp.plot_detector(t_trg_met,geometry)
	smp.plot_galactic_plane()
	smp.plot_sum(t_trg_met,geometry)
	smp.plot_moon(t_trg_met,geometry)
	smp.plot_continue_source()
	smp.savefig(savedir_i + 'A_sky_map.png')
	smp.close()

	if location50_300[0] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location50_300[0]
		c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savedir_i + 'B_location_50_300_chi2.csv',index=False)

		c1 = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
		save_perspec1 = pd.DataFrame(c1)
		save_perspec1.to_csv(savedir_i +'B_location_50_300.csv',index=False)

	if location5_50[0] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location5_50[0]
		c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savedir_i + 'B_location_5_50_chi2.csv',index=False)
		c1 = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
		save_perspec1 = pd.DataFrame(c1)
		save_perspec1.to_csv(savedir_i +'B_location_5_50.csv',index=False)

	lc_arri = lc_arr[0].T
	bs_arri = bs_arr[0].T
	# SNR_arri = SNR_arr[2].T
	lc_t = bins_lightcurve_c
	fig = plt.figure(figsize=(10, 10))
	fig.suptitle('lc')
	for j in range(len(lc_arri)):
		detec = name[j]
		plt.subplot(6, 2, j + 1)
		plt.plot(lc_t - t_trg_met, lc_arri[j], label='lc ' + detec)
		plt.plot(lc_t - t_trg_met, bs_arri[j], label='bs ' + detec)

		plt.axvline(x=t1_met - t_trg_met, color='r')
		plt.axvline(x=t2_met - t_trg_met, color='g')
		plt.legend()
	plt.savefig(savedir_i + 'A_lc.png')
	plt.close()

	lc_arri = lc_arr[1].T
	bs_arri = bs_arr[1].T
	# SNR_arri = SNR_arr[2].T
	lc_t = bins_lightcurve_c
	fig = plt.figure(figsize=(10, 10))
	fig.suptitle('lc 5-50kev')
	for j in range(len(lc_arri)):
		detec = name[j]
		plt.subplot(6, 2, j + 1)
		plt.plot(lc_t - t_trg_met, lc_arri[j], label='lc ' + detec)
		plt.plot(lc_t - t_trg_met, bs_arri[j], label='bs ' + detec)
		plt.axvline(x=t1_met - t_trg_met, color='r')
		plt.axvline(x=t2_met - t_trg_met, color='g')
		plt.legend()
	plt.savefig(savedir_i + 'A_lc_5_50.png')
	plt.close()

	lc_arri = lc_arr[2].T
	bs_arri = bs_arr[2].T
	# SNR_arri = SNR_arr[2].T
	lc_t = bins_lightcurve_c
	fig = plt.figure(figsize=(10, 10))
	fig.suptitle('lc 50-300kev')
	for j in range(len(lc_arri)):
		detec = name[j]
		plt.subplot(6, 2, j + 1)
		plt.plot(lc_t - t_trg_met, lc_arri[j], label='lc ' + detec)
		plt.plot(lc_t - t_trg_met, bs_arri[j], label='bs ' + detec)
		plt.axvline(x=t1_met - t_trg_met, color='r')
		plt.axvline(x=t2_met - t_trg_met, color='g')
		plt.legend()
	plt.savefig(savedir_i + 'A_lc_50_300.png')
	plt.close()





