import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from astropy.stats import sigma_clip, mad_std
import operator
from scipy import stats


savedir = '/home/laojin/my_lat/DV/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

rate = 15000
time_start = 0
time_stop = 600
binsize1 = 0.01
binsize2 = 1.0

sample_start = -10
sample_stop = 610

np.random.seed(115)#固定随机数
T = sample_stop-sample_start
#t = stats.uniform.rvs(size = int(rate*T))*T-10
t = np.random.uniform(sample_start,sample_stop,rate*T)
y = np.random.uniform(0,10,rate*T)

edge1 = np.arange(time_start,time_stop,binsize1)
edge2 = np.arange(time_start,time_stop,binsize2)

bin_n1, edges1 = np.histogram(t,bins = edge1)
t_c1 = (edges1[1:]+edges1[:-1])*0.5
rate1 = bin_n1/binsize1
mean1 = rate1.mean()
sigma1 = rate1.std(ddof = 1)#必须考虑有偏

mask1 = sigma_clip(rate1,sigma=5,maxiters=5,stdfunc=mad_std).mask
myfilter = list(map(operator.not_, mask1))
rate1_median_part = rate1[myfilter]
part_bins1 = np.linspace(rate1.min(),rate1.max(),100)
histvalue1,histbin1 = np.histogram(rate1,bins = part_bins1)
histvalue1 = np.concatenate((histvalue1,histvalue1[-1:]))
loc1,scale1 = stats.norm.fit(rate1_median_part)
Y1 = stats.norm(loc=loc1,scale=scale1)



bin_n2, edges2 = np.histogram(t,bins = edge2)
t_c2 = (edges2[1:]+edges2[:-1])*0.5
rate2 = bin_n2/binsize2
mean2 = rate2.mean()
sigma2 = rate2.std(ddof = 1)#必须考虑有偏

mask2 = sigma_clip(rate2,sigma=5,maxiters=5,stdfunc=mad_std).mask
myfilter2 = list(map(operator.not_, mask2))
rate2_median_part = rate2[myfilter2]
part_bins2 = np.linspace(rate1.min(),rate1.max(),100)
histvalue2,histbin2 = np.histogram(rate2,bins = part_bins2)
histvalue2 = np.concatenate((histvalue2,histvalue2[-1:]))
loc2,scale2 = stats.norm.fit(rate2_median_part)
Y2 = stats.norm(loc=loc2,scale=scale2)

fig = plt.figure(figsize = (10,10))
gs = mpl.gridspec.GridSpec(3, 4,wspace=0.05, hspace=0.05)
ax0 = fig.add_subplot(gs[0,:3])
ax0.plot(t,y,',',markersize = 0.001,color = 'k')
ax0.set_ylabel('A part of sample space')
ax0.set_xlim([0,600])
ax0.set_ylim([0,0.05])
ax0.set_xticks([])
ax0.set_yticks([])

ax1 = fig.add_subplot(gs[1,:3])
ax1.plot(t_c1,rate1,label = r'Rate curve of dt = %.2fs' % binsize1)
ax1.set_ylim([10000,20000])
ax1.set_xlim([0,600])
ax1.set_xticks([])
ax1.axhline(y = mean1,color = 'r',label = r'$Mean_{dt = %.2fs} = %.2f$' % (binsize1,mean1))
ax1.axhline(y = sigma1+mean1,color = '#f47920',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize1,sigma1))#大
ax1.axhline(y = mean1-sigma1,color = '#f47920')
ax1.set_ylabel('rate')
ax1.legend()


ax2 = fig.add_subplot(gs[2,:3])
ax2.plot(t_c2,rate2,label = r'Rate curve of dt = %.2fs' % binsize2)
ax2.set_ylim([10000,20000])
ax2.set_xlim([0,600])
ax2.axhline(y = mean2,color = 'r',label = r'$Mean_{dt = %.2fs} = %.2f$' % (binsize2,mean2))
ax2.axhline(y = sigma2+mean2,color = 'g',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize2,sigma2))
ax2.axhline(y = mean2-sigma2,color = 'g')
ax2.axhline(y = mean1+sigma1,color = '#f47920',label =r'$\sigma_{dt = %.2fs}$ = %.2f' % (binsize1,sigma1))
ax2.axhline(y = mean1-sigma1,color = '#f47920')
ax2.set_xlabel('t')
ax2.set_ylabel('rate')
ax2.legend()

#ax3 = fig.add_subplot(gs[1,2], sharey=ax1)
ax3 = fig.add_subplot(gs[1,3])
ax3.fill_betweenx(histbin1,histvalue1,0,step= 'post')
ax3.plot(Y1.pdf(part_bins1)*rate1.size*(part_bins1[1]-part_bins1[0]),part_bins1,color = '#f47920')
ax3.set_xlim(0)
ax3.set_ylim([10000,20000])
ax3.set_xticks([])
ax3.set_yticks([])

ax4 = fig.add_subplot(gs[2,3])
ax4.fill_betweenx(histbin2,histvalue2,0,step= 'post')
ax4.plot(Y2.pdf(part_bins2)*rate2.size*(part_bins2[1]-part_bins2[0]),part_bins2,color = '#f47920')
ax4.set_xlim(0)
ax4.set_ylim([10000,20000])
ax4.set_xticks([])
ax4.set_yticks([])

fig.savefig(savedir + 'A_sample_rate1.png')
plt.close(fig)


'''
plt.figure(figsize = (10,10))
plt.subplot(3,1,1)
plt.plot(t,y,',',markersize = 0.001,color = 'k')
plt.ylabel('A part of sample space')
plt.yticks([])
plt.xlim([0,600])
plt.ylim([0,0.05])
plt.subplot(3,1,2)
plt.plot(t_c1,rate1,label = r'rate curve of dt = %.2fs' % binsize1)
plt.ylim([10000,20000])
plt.xlim([0,600])
plt.axhline(y = mean1,color = 'r',label = r'$mean_{dt = %.2fs} = %.2f$' % (binsize1,mean1))
plt.axhline(y = sigma1+mean1,color = '#f47920',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize1,sigma1))#大
plt.axhline(y = mean1-sigma1,color = '#f47920')
plt.ylabel('rate')
plt.legend()

plt.subplot(3,1,3)
plt.plot(t_c2,rate2,label = r'rate curve of dt = %.2fs' % binsize2)
plt.ylim([10000,20000])
plt.xlim([0,600])
plt.axhline(y = mean2,color = 'r',label = r'$mean_{dt = %.2fs} = %.2f$' % (binsize2,mean2))
plt.axhline(y = sigma2+mean2,color = 'g',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize2,sigma2))
plt.axhline(y = mean2-sigma2,color = 'g')
plt.axhline(y = mean1+sigma1,color = '#f47920',label =r'$\sigma_{dt = %.2fs}$ = %.2f' % (binsize1,sigma1))
plt.axhline(y = mean1-sigma1,color = '#f47920')
plt.xlabel('t')
plt.ylabel('rate')
plt.legend()
plt.savefig(savedir + 'A_sample_rate1.png')
plt.close()
'''
plt.figure()
plt.plot(t_c2,rate2,label = r'Rate curve of dt = %.2fs' % binsize2)
plt.xlim([0,600])
plt.axhline(y = mean2,color = 'r',label = r'$Mean_{dt = %.2fs} = %.2f$' % (binsize2,mean2))
plt.axhline(y = sigma2+mean2,color = 'g',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize2,sigma2))
plt.axhline(y = mean2-sigma2,color = 'g')
sigma2_1 = np.sqrt(binsize2)*sigma2/np.sqrt(binsize1)
plt.axhline(y = mean2+sigma2_1,color = 'k',label = r'$\sigma_{%.2f \Rightarrow %.2fs}$ = %.2f' % (binsize2,binsize1,sigma2_1))
plt.axhline(y = mean1+sigma1,color = '#f47920',label = r'$\sigma_{dt = %.2fs}$ = %.2f' % (binsize1,sigma1))
plt.axhline(y = mean2-sigma2_1,color = 'k')
plt.axhline(y = mean1-sigma1,color = '#f47920')
plt.xlabel('t')
plt.ylabel('rate')
plt.legend()
plt.savefig(savedir + 'A_sample_rate2.png')
plt.close()







