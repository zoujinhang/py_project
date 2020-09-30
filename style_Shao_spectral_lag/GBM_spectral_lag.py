
import numpy as np
import matplotlib.pyplot as plt
from .tool import TD_bs,WhittakerSmooth,findfile,ch_to_energy,get_band
from scipy.interpolate import interp1d
import os
from astropy.io import fits
from matplotlib.gridspec import GridSpec
import pandas as pd




class GRB(object):
	
	def __init__(self,bnname,database,savedir):
		
		year = '20'+bnname[2:4]
		self.NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
		self.BGO = ['b0','b1']
		self.link = database + year + '/' + bnname + '/'
		self.savetop = savedir
		self.bnname = bnname
		
	def get_burst_data(self):
		
		data = {}
		for ni in self.NaI:
			ni_list = findfile(self.link, 'glg_tte_' + ni + '_' + self.bnname+ '_v*')
			if len(ni_list) >= 1:
				# print(ni)
				data[ni] = fits.open(self.link + ni_list[0])
			else:
				data[ni] = None
		for bi in self.BGO:
			bi_list = findfile(self.link, 'glg_tte_' + bi + '_' + self.bnname + '_v*')
			if len(bi_list) >= 1:
				data[bi] = fits.open(self.link + bi_list[0])
			else:
				data[bi] = None
				
		return 	data
	
	def get_lag(self,detector,edges,wind,binsize = 0.2,sigma = 4,Positive_lag_positive_value = True):
		
		band = get_band([10,1000],10,ovelap=0.5)
		plot_savedir = self.savetop + 'A_check/'
		
		files = self.get_burst_data()
		hl = files[detector]
		trigtime = hl[0].header['TRIGTIME']
		time = hl[2].data.field(0)
		ch = hl[2].data.field(1)
		ch_n = hl[1].data.field(0)
		e1 = hl[1].data.field(1)
		e2 = hl[1].data.field(2)
		t = time - trigtime
		
		bins = np.arange(min(edges),max(edges),binsize)
		t,energy = ch_to_energy(t,ch,ch_n,e1,e2)
		
		t_c = 0.5 * (bins[1:] + bins[:-1])
		if wind is not None:
			wind_index = np.where((t_c >= wind[0]) & (t_c <= wind[-1]))[0]
		else:
			wind_index = np.arange(0, len(t_c), dtype=int)
		event_band = get_event_with_band(t, energy, band)
		cs0 = None
		cs_err0 = None
		lag_errh2 = 0
		lag_errl2 = 0
		lag_all = 0
		return_list = []
		lc_list = []
		for index_, (t_i, e_i) in enumerate(event_band):
			if cs0 is None or cs_err0 is None:
				num_i = np.histogram(t_i, bins=bins)[0]
				num_err_i = np.sqrt(num_i)
				rate_i = num_i / binsize
				rate_err_i = num_err_i / binsize
				csi, bsi, sigma_i = TD_bs(t_c, rate_i, sigma=True)
				lc_list.append([rate_i, bsi, sigma_i])
				SNRi = csi / sigma_i
				if SNRi.max() > sigma:
					cs0 = csi
					cs_err0 = rate_err_i
					return_list.append([index_, 0, 1, 1])
			else:
				num_i = np.histogram(t_i, bins=bins)[0]
				num_err_i = np.sqrt(num_i)
				rate_i = num_i / binsize
				rate_err_i = num_err_i / binsize
				csi, bsi, sigma_i = TD_bs(t_c, rate_i, sigma=True)
				lc_list.append([rate_i, bsi, sigma_i])
				SNRi = csi / sigma_i
				if SNRi.max() > sigma:
					if plot_savedir is not None:
						if os.path.exists(plot_savedir) == False:
							os.makedirs(plot_savedir)
						savename = plot_savedir + 'A_nccf_' + str(index_) + '.png'
					else:
						savename = None
					lag, lag_errl, lag_errh = get_one_lag(cs0[wind_index], csi[wind_index],
					                                      cs_err0[wind_index], rate_err_i[wind_index],
					                                      t_c[wind_index], mcmc_num=500, save=savename)
					lag_all = lag_all + lag
					lag_errl2 = lag_errl2 + lag_errl ** 2
					lag_errh2 = lag_errh2 + lag_errh ** 2
					return_list.append([index_, lag_all, np.sqrt(lag_errl2), np.sqrt(lag_errh2)])
					cs0 = csi
					cs_err0 = rate_err_i
		lag0 = np.array(return_list).T
		if Positive_lag_positive_value:
			index_,lag,lag_errl,lag_errh = lag0
		else:
			index_,lag,lag_errl,lag_errh = np.vstack([lag0[0],-1*lag0[1],lag0[3],lag0[2]])
		n = len(lc_list)
		fig = plt.figure(constrained_layout=True,figsize = (5,2*n))
		gs = GridSpec(n, 1, figure=fig)
		for index_i,(ratei,bsi,sigmai) in enumerate(lc_list):
			bandi = band[index_i]
			labeli = '%.1f-%.1f keV' % (bandi[0],bandi[-1])
			ax = fig.add_subplot(gs[-index_i-1])
			ax.plot(t_c,ratei,'-',color = 'k',label = labeli)
			ax.plot(t_c,bsi,'-',color = 'r')
			ax.plot(t_c,bsi+sigma*sigmai,'-',color = 'g')
			if wind is not None:
				ax.set_xlim(wind[0],wind[-1])
			else:
				ax.set_xlim(t_c.min(),t_c.max())
			ax.set_ylabel('Rate')
			ax.legend()
			if index_i != 0 :
				ax.tick_params(labelbottom = False)
			else:
				ax.set_xlabel('Time')
		plt.savefig(self.savetop+'B_lightcurve.png')
		plt.close()
		band_l, band_h = (band[index_.astype(np.int)]).T
		energyc = np.sqrt(band_l*band_h)
		
		fig = plt.figure(constrained_layout=True)
		ax = fig.add_subplot(1,1,1)
		ax.errorbar(energyc[1:],lag[1:],yerr = [lag_errl[1:],lag_errh[1:]],elinewidth=1,capsize=2,label = 'lag data',fmt = '.')
		ax.set_ylabel('lag (s)')
		ax.set_xlabel('energy Kev')
		plt.savefig(self.savetop+'C_lag.png')
		plt.close()
		
		c = {'band_l':band_l,
		     'band_h':band_h,
		     'lag':lag,
		     'lag_errl':lag_errl,
		     'lag_errh':lag_errh}
		hl = pd.DataFrame(c)
		hl.to_csv(self.savetop+'D_lag.csv',index=False)



def get_one_lag(cs1,cs2,cs1_err,cs2_err,time,mcmc_num = 1000,save = None):
	
	n = len(time)
	dt = np.abs(time[1]-time[0])
	lag_t = np.arange(-n+1,n,1)*dt
	nccf = np.correlate(cs1 / cs1.max(), cs2 / cs2.max(), 'full')
	lag_time = get_nccf_peak_time(lag_t,nccf,save = save)
	lag_list = []
	for i in range(mcmc_num):
		rand1 = np.random.randn(n)
		cs1_sim = cs1 + cs1_err*rand1
		rand2 = np.random.randn(n)
		cs2_sim = cs2 + cs2_err*rand2
		nccf_sim = np.correlate(cs1_sim/cs1_sim.max(),cs2_sim/cs2_sim.max(),'full')
		lag_sim = get_nccf_peak_time(lag_t,nccf_sim)
		lag_list.append(lag_sim)
	lag_list = np.array(lag_list)
	lag_mean = lag_list.mean()
	lag_sigma = lag_list.std(ddof = 1)
	dx = lag_mean-lag_time
	lag_errl = lag_sigma-dx
	lag_errh = lag_sigma+dx
	
	return lag_time,lag_errl,lag_errh
	
	
def get_nccf_peak_time(lag_t,nccf,precision = 0.001,save = None):
	
	dt = np.abs(lag_t[1]-lag_t[0])
	w = np.ones(len(nccf))
	smoo_nccf = WhittakerSmooth(nccf,w,0.15/dt**1.5)
	intef = interp1d(lag_t,smoo_nccf,kind = 'quadratic')
	new_lat_t = np.arange(lag_t[1],lag_t[-2]+precision,precision)
	new_nccf = intef(new_lat_t)
	index_ = np.argmax(new_nccf)
	t_max = new_lat_t[index_]
	if save is not None:
		fig = plt.figure(constrained_layout=True)
		ax = fig.add_subplot(1,1,1)
		ax.plot(lag_t,nccf,'.',color = 'k',label = 'nccf')
		ax.plot(new_lat_t,new_nccf,'-',color='#f47920',label = 'fit')
		ax.axvline(x =t_max,color = 'r' )
		ax.set_xlim(t_max-5,t_max+5)
		ax.set_xlabel('lag time (s)')
		ax.set_ylabel('nccf')
		ax.legend()
		fig.savefig(save)
		plt.close(fig)
	return t_max

def get_event_with_band(time,energy,band):
	
	return_list = []
	
	for (el,eh) in band:
		index = np.where((energy>=el)&(energy<=eh))[0]
		t_i = time[index]
		e_i = energy[index]
		return_list.append([t_i,e_i])
	
	return return_list



















