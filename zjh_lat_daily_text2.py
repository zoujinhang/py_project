import numpy as np
import matplotlib.pyplot as plt
from Fermi_tool.daily import Database
import os
from Data_analysis import TD_baseline,background_correction,get_bayesian_duration,WhittakerSmooth
from astropy.stats import bayesian_blocks,sigma_clip,mad_std
import operator
from scipy import stats
import warnings


savedir = '/home/laojin/my_lat/daily_text/check/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)

def get_w(cs,sigma):
	return np.exp(-0.5/(sigma**2)*(cs)**2)

def TD_bs(t,rate,it = 2,lambda_=4000,sigma = False):
	dt = t[1]-t[0]
	t_c,cs,bs = TD_baseline(t,rate)
	scale = 0
	for i in range(it):
		mask = sigma_clip(cs,sigma=5,maxiters=5,stdfunc=mad_std).mask
		myfilter = list(map(operator.not_, mask))
		lc_median_part = cs[myfilter]
		loc,scale = stats.norm.fit(lc_median_part)
		w = get_w(cs,scale)
		bs = WhittakerSmooth(rate,w,lambda_=lambda_/dt**1.5)
		cs = rate - bs
	if sigma:
		return cs,bs,scale
	else:
		return cs,bs

def zero_check(t,r):
	pass


topdir = '/media/laojin/Elements/daily/'

binsize = 0.064
#t_start = '2020-02-19T11:35:22.594'
#t_stop = '2020-02-19T11:35:43.546'
t_start = '2013-04-27T07:46:30'
t_stop = '2013-04-27T07:48:00'
fermi = Database(topdir)

data = fermi.get_detector_data(t_start,t_stop)

ni = data['n6']['events']
time = ni['TIME'].values

t0 = time.min()
t = time - t0

bins = np.arange(t.min(),t.max(),binsize)
t_c = 0.5*(bins[:-1]+bins[1:])
bin_n = np.histogram(t,bins = bins)[0]
rate = bin_n/binsize

t_c,cs,bs = TD_baseline(t_c,rate)
cs1,bs1,sigma =TD_bs(t_c,rate,sigma = True)


new_rate = cs + bs.mean()
warnings.simplefilter('error')
try:
	n_index = np.where(new_rate>0)[0]
	nn_new_rate = new_rate[n_index]
	nn_t_c = t_c[n_index]
	edges = bayesian_blocks(nn_t_c,np.round(nn_new_rate*binsize),fitness = 'events',gamma = np.exp(-5))
	results = background_correction(t_c, new_rate, edges, degree=5, plot_save=savedir + 'A_.png')
	startedges, stopedges, new_snr = get_bayesian_duration(results, sigma=4, max_snr=True)
	print(startedges, stopedges, new_snr)
	t_c1, new_rate1 = results['lc']
	re_rate = results['re_hist'][0]
	re_rate = np.concatenate((re_rate[:1], re_rate))
	
	plt.plot(t_c1, new_rate1)
	plt.step(edges, re_rate)
	plt.savefig(savedir + 'B_.png')
	plt.close()
	
	plt.plot(t_c, rate)
	plt.plot(t_c, bs)
	plt.plot(t_c, bs + 3 * sigma)
	plt.plot(t_c, bs1, color='r')
	plt.plot(t_c, bs1 + 3 * sigma, color='r')
	plt.savefig(savedir + 'C_.png')
	plt.close()
	
except(RuntimeWarning):
	
	print('encountered')
warnings.simplefilter('always')



