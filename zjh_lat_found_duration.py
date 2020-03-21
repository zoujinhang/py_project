import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from astropy.stats import bayesian_blocks
from Calculate_duration.background_kernel import r_baseline

datalink = '/home/laojin/trigdata/glg_tte_n1_bn200101861_v00.fit'
savedir = '/home/laojin/my_lat/found_duration/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

def re_histogram(t,rate,edges):
	'''
	now it is a stable function, it can re-bin the data, 't' and 'rate', with the new bin edges, 'edges'.
	:param t: The time of light curve.
	:param rate: The counts rate of light curve.
	:param edges: the new bin-edges, whose interval can different from each other.
	:return: two array ,The rebined rate and sigma of each bins.
	'''
	index = np.where(t>=edges[0])[0]
	t = t[index]
	rate = rate[index]
	re_rate = np.zeros(len(edges)-1)
	re_sigma = np.zeros(len(edges)-1)
	edges_num = 1
	index_list = []
	ones = []
	for index1,value in enumerate(t):
		if value > edges[edges_num] or value == t[-1]:
			index_list.append(ones)
			ones = []
			edges_num = edges_num+1
			if edges_num == len(edges):
				break
		else:
			ones.append(index1)
	for index1,value in enumerate(index_list):
		one_rate = rate[value]
		re_rate[index1] = one_rate.mean()
		if len(one_rate) == 1:
			re_sigma[index1] = one_rate.std(ddof = 0)
		else:
			re_sigma[index1] = one_rate.std(ddof = 1)
	return re_rate,re_sigma





hl = fits.open(datalink)

trigtime = hl[0].header['TRIGTIME']

time = hl[2].data.field(0)
ch = hl[2].data.field(1)

t = time - trigtime

binsize = 0.064
bins = np.arange(t[0],t[-1],binsize)
bin_n,bin_edges = np.histogram(t,bins = bins)

t_c = (bin_edges[1:] + bin_edges[:-1])*0.5
rate = bin_n/binsize
bs_rate = r_baseline(rate,binsize)

rate_sm = rate - bs_rate+bs_rate.mean()
bin_n_sm = np.round(rate_sm*binsize)

edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',p0 = 0.001)
print('edges:\n',edges)


re_rate,re_sigma = re_histogram(t_c,rate_sm,edges)
re_rate1 = np.concatenate((re_rate[:1],re_rate))
re_sigma1 = np.concatenate((re_sigma[:1],re_sigma))

plt.figure()
plt.subplot(2,1,1)
plt.plot(t_c,rate)
plt.plot(t_c,bs_rate)
plt.ylim(0,5000)
plt.xlim(edges[0],edges[-1])
plt.subplot(2,1,2)
plt.plot(t_c,rate_sm)
for i in edges:
	plt.axvline(x = i,color = 'r')
plt.step(edges,re_rate1,color = 'k')
plt.step(edges,re_rate1+re_sigma1,color = 'y')
plt.step(edges,re_rate1-re_sigma1,color = 'y')

plt.ylim(0,5000)
plt.xlim(edges[0],edges[-1])
plt.savefig(savedir + 'A_lightcurve.png')
plt.close()










