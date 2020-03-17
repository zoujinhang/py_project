import numpy as np
import matplotlib.pyplot as plt
import os

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

T = sample_stop-sample_start
t = np.random.uniform(sample_start,sample_stop,rate*T)
y = np.random.uniform(0,10,rate*T)

edge1 = np.arange(time_start,time_stop,binsize1)
edge2 = np.arange(time_start,time_stop,binsize2)

bin_n1, edges1 = np.histogram(t,bins = edge1)
t_c1 = (edges1[1:]+edges1[:-1])*0.5
rate1 = bin_n1/binsize1
mean1 = rate1.mean()
sigma1 = rate1.std(ddof = 1)#必须考虑有偏

bin_n2, edges2 = np.histogram(t,bins = edge2)
t_c2 = (edges2[1:]+edges2[:-1])*0.5
rate2 = bin_n2/binsize2
mean2 = rate2.mean()
sigma2 = rate2.std(ddof = 1)#必须考虑有偏

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

plt.figure()
plt.plot(t_c2,rate2,label = r'rate curve of dt = %.2fs' % binsize2)
plt.xlim([0,600])
plt.axhline(y = mean2,color = 'r',label = r'$mean_{dt = %.2fs} = %.2f$' % (binsize2,mean2))
plt.axhline(y = sigma2+mean2,color = 'g',label = r'$\sigma_{dt = %.2fs} = %.2f$' % (binsize2,sigma2))
plt.axhline(y = mean2-sigma2,color = 'g')
sigma2_1 = np.sqrt(binsize2)*sigma2/np.sqrt(binsize1)
plt.axhline(y = mean2+sigma2_1,color = 'k',label = r'$\sigma_{dt = %.2f \Rightarrow %.2fs}$ = %.2f' % (binsize2,binsize1,sigma2_1))
plt.axhline(y = mean1+sigma1,color = '#f47920',label = r'$\sigma_{dt = %.2fs}$ = %.2f' % (binsize1,sigma1))
plt.axhline(y = mean2-sigma2_1,color = 'k')
plt.axhline(y = mean1-sigma1,color = '#f47920')
plt.xlabel('t')
plt.ylabel('rate')
plt.legend()
plt.savefig(savedir + 'A_sample_rate2.png')
plt.close()







