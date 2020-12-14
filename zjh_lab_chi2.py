import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os

savedir = '/home/laojin/my_lat/chi2_distribution/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

x = np.linspace(0,60,100)
y = stats.chi2.pdf(x,df=12)

plt.plot(x,y)
plt.savefig(savedir + 'A_chi2_12.png')
plt.close()


xm = np.random.randn(10000)
ym = np.random.randn(10000)*2
zm = np.random.randn(10000)*2

r2 = np.zeros_like(xm)
q = np.zeros_like(xm)
for i in range(xm.size):

	r2[i] = ((xm-xm[i])**2+(ym-ym[i])**2+(zm-zm[i])**2).sum()
	q[i] = (1/(np.sqrt((xm-xm[i])**2+(ym-ym[i])**2+(zm-zm[i])**2)+1)).sum()
r2 = np.sqrt(r2)

#q = 1/(r2+1)
r2max = r2.max()
r2min = r2.min()
stad = np.linspace(r2min,r2max,20)
plt.figure(figsize = (5,5))
plt.plot(xm,ym,',',color = 'k',markersize = 2)
#plt.tricontourf(xm,ym,r2,stad)
plt.tricontour(xm,ym,r2,stad)
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.savefig(savedir + 'A_dis.png')
plt.close()

r2max = q.max()
r2min = q.min()
stad = np.linspace(r2min,r2max,20)
plt.figure(figsize = (5,5))
plt.plot(xm,ym,',',color = 'k',markersize = 2)
#plt.tricontourf(xm,ym,q,stad)
plt.tricontour(xm,ym,q,stad)
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.savefig(savedir + 'A_q.png')
plt.close()
