import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from astropy.stats import bayesian_blocks
from Data_analysis import TD_baseline,background_correction,get_bayesian_duration,get_SNR


datalink = '/home/laojin/trigdata/glg_tte_n1_bn200101861_v00.fit'
savedir = '/home/laojin/my_lat/bb_duration_text/'

if os.path.exists(savedir) == False:
	os.makedirs(savedir)

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
t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)

rate_sm = cs_rate+bs_rate.mean()
bin_n_sm = np.round(rate_sm*binsize)

edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',p0 = 0.001)
#-----------------------------------------------------------------
result = background_correction(t_c,rate_sm,edges,degree = 50)
startedges,stopedges = get_bayesian_duration(result,sigma = 5)
#-----------------------------------------------------------------
background_mean,background_sigma,backgound_size = result['bkg']
t_c,rate_correct = result['lc']
re_rate,re_sigma,index_list = result['re_hist']

re_rate1 = np.concatenate((re_rate[:1],re_rate))
re_sigma1 = np.concatenate((re_sigma[:1],re_sigma))

SNR = get_SNR(edges,re_rate,background_mean,background_sigma,binsize)
SNR1 = np.concatenate((SNR[:1],SNR))

plt.figure()
plt.subplot(2,1,1)
plt.plot(t_c,rate)
plt.plot(t_c,bs_rate)
plt.ylim(0,5000)
plt.xlim(edges[0],edges[-1])
plt.subplot(2,1,2)
plt.plot(t_c,rate_sm)

plt.step(edges,re_rate1,color = 'k')
#plt.step(edges,re_rate1+re_sigma1,color = 'y')
#plt.step(edges,re_rate1-re_sigma1,color = 'y')
plt.axhline(y = background_mean,color = 'r')
#plt.axhline(y = background_mean-background_sigma,color = 'r')
#plt.axhline(y = background_mean+background_sigma,color = 'r')
plt.ylim(0,5000)
plt.xlim(edges[0],edges[-1])
plt.savefig(savedir + 'A_lightcurve.png')
plt.close()

fig = plt.figure()
gs = mpl.gridspec.GridSpec(2,1,wspace = 0.05,hspace = 0.05)
ax1 = fig.add_subplot(gs[0,0])
ax1.plot(t_c,rate_correct,label = 'lightcurve')
ax1.step(edges,re_rate1,color = 'k',label = 'Bayesian blocks')
for v in startedges:
	ax1.axvline(x = v,color = 'r')
for v in stopedges:
	ax1.axvline(x = v,color = 'g')
ax1.set_xticks([])
ax1.set_ylabel('Rate')
ax1.set_xlim(-100,200)
ax1.set_ylim(0,5000)
ax1.legend()

ax2 = fig.add_subplot(gs[1,0])
ax2.plot(t_c,(rate_correct-background_mean)/background_sigma,label = 'Classical SNR')
ax2.step(edges,SNR1,color = 'k',label = 'Bayesian blocks SNR')
for v in startedges:
	ax2.axvline(x = v,color = 'r')
for v in stopedges:
	ax2.axvline(x = v,color = 'g')
ax2.axhline(y=5,color = 'r',label = r'5 $\sigma$')
ax2.axhline(y=0,color = 'y',label = 'Zero')
ax2.set_xlabel('time (s)')
ax2.set_ylabel('SNR')
ax2.set_xlim(-100,200)
ax2.set_ylim(-10,100)
ax2.legend()

fig.savefig(savedir + 'C_SNR_lightcurve.png')
plt.close(fig)








