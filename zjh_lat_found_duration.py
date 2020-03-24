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
	
	index = np.where(t >= edges[0])[0]
	t = t[index]
	rate = rate[index]
	re_rate = np.zeros(len(edges)-1)
	re_sigma = np.zeros(len(edges)-1)
	edges_num = 1
	index_list = []
	ones = []
	for index1,value in enumerate(t):
		if value > edges[edges_num] or value == t[-1]:
			ones.append(index1)
			index_list.append(ones)
			ones = []
			edges_num = edges_num+1
			if edges_num == len(edges):
				break
		else:
			ones.append(index1)
	for index1,value in enumerate(index_list):
		one_rate = rate[value]
		index_list[index1] = np.array(value)+index[0]
		re_rate[index1] = one_rate.mean()
		if len(one_rate) == 1:
			re_sigma[index1] = one_rate.std(ddof = 0)
		else:
			re_sigma[index1] = one_rate.std(ddof = 1)
	return re_rate,re_sigma,np.array(index_list)

def confidence_analysis(mean1,sigma1,T1,mean2,sigma2,T2,dt,degree = 3):
	
	'''
	conditional dependence of confidence
	:param mean1:
	:param sigma1:
	:param T1:
	:param mean2:
	:param sigma2:
	:param T2:
	:param dt:
	:param degree:
	:return:
	'''
	
	sigma1_2 = sigma1*np.sqrt(dt/T2)
	sigma2_1 = sigma2*np.sqrt(dt/T1)
	return (mean1>=mean2-degree*sigma2_1)&(mean1<=mean2+degree*sigma2_1)&(mean2>=mean1-degree*sigma1_2)&(mean2<=mean1+degree*sigma1_2)


def found_background(t,rate,index_list,edges,re_rate = None ,sigma = None,degree = 70):
	'''
	
	:param t: Time of light curve.
	:param rate: counts rate of light curve.
	:param index_list: the index of light curve in each block bettween each edge.
	:param edges: edges obtained from Bayesian bloacks
	:param re_rate: the rate in each block.
	:param sigma: the sigma in each block.
	:param degree: Tolerance of background consistency. Default value is 70.
	:return: three value of mean, sigma, and size of background.
	'''
	
	if re_rate is None or sigma is None:
		sigma = np.zeros(len(edges)-1)
		re_rate = np.zeros(len(edges)-1)
		for index1,value in enumerate(index_list):
			one_rate = rate[value]
			re_rate[index1] = one_rate.mean()
			if len(one_rate) == 1:
				sigma[index1] = one_rate.std(ddof = 0)
			else:
				sigma[index1] = one_rate.std(ddof = 1)
	dt = t[1]-t[0]
	binsize = edges[1:]-edges[:-1]
	sort_index = np.argsort(-binsize)
	sort_binsize = binsize[sort_index]
	sort_sigma = sigma[sort_index]
	sort_re_rate = re_rate[sort_index]
	sort_index_list = np.array(index_list)[sort_index]
	
	background_pool = rate[sort_index_list[0]]
	mean1 = sort_re_rate[0]
	sigma1 = sort_sigma[0]
	binsize1 = sort_binsize[0]
	n = 0
	for i in range(1,len(sort_binsize)):
		
		if confidence_analysis(mean1,sigma1,binsize1,sort_re_rate[i],sort_sigma[i],sort_binsize[i],dt,degree):
			n = n+1
			background_pool = np.concatenate((background_pool,rate[sort_index_list[i]]))
			mean1 = background_pool.mean()
			sigma1 = background_pool.std(ddof = 1)
			binsize1 = binsize1 + sort_binsize[i]
	print('The number of background blocks is %d.'%n)
	return mean1,sigma1,binsize1



def core_of_get_SNR(edges,re_rate,background_mean,background_sigma,dt):
	binsize = edges[1:] - edges[:-1]
	sigma = background_sigma * np.sqrt(dt/binsize)
	SNR = (re_rate - background_mean)/sigma
	SNR[SNR<0] = 0 #we only care about the part where the SNR is greater than 0.
	return SNR



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


re_rate,re_sigma,index_list = re_histogram(t_c,rate_sm,edges)

background_mean,background_sigma,backgound_size = found_background(t_c,rate_sm,index_list,edges,re_rate,re_sigma,degree = 70)

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

plt.step(edges,re_rate1,color = 'k')
plt.step(edges,re_rate1+re_sigma1,color = 'y')
plt.step(edges,re_rate1-re_sigma1,color = 'y')
plt.axhline(y = background_mean,color = 'r')
plt.axhline(y = background_mean-background_sigma,color = 'r')
plt.axhline(y = background_mean+background_sigma,color = 'r')
plt.ylim(0,5000)
plt.xlim(edges[0],edges[-1])
plt.savefig(savedir + 'A_lightcurve.png')
plt.close()

SNR = core_of_get_SNR(edges,re_rate,background_mean,background_sigma,binsize)
SNR1 = np.concatenate((SNR[:1],SNR))
'''
the classical SNR is different from bayesian SNR . because bayesian SNR has different value of sigma.
'''

plt.figure()
plt.plot(t_c,(rate_sm-background_mean)/background_sigma,label = 'classical SNR')
plt.step(edges,SNR1,color = 'k',label = 'Bayesian blocks SNR')
plt.axhline(y=5,color = 'r',label = r'5 $\sigma$')
plt.axhline(y=0,color = 'y',label = 'zero')
plt.xlim(edges[0],edges[-1])
plt.xlabel('time (s)')
plt.ylabel('SNR')
plt.legend()
plt.savefig(savedir+ 'B_SNR.png')
plt.close()






