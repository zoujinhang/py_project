import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import os


def K(s):
	return mp.besselk(5/3,s)

def F(x):
	f = np.zeros(len(x))
	for index,value in enumerate(x):
		f[index] = value*mp.quad(K,[value,np.inf])
	return f

savedir = '/home/laojin/my_lat/synchrotron_radiation/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

#a  = mp.besselk(5/3,0)
#print(a)

x = np.logspace(-3,1,100)
y = F(x)

plt.plot(x,y,color = 'k',label = 'F(x)')
plt.axvline(x = 0.3,color = 'r',label = 'x = 0.3')
plt.ylabel('F(x)')
plt.xlabel('x')
plt.xscale('log')
plt.xlim(x[0],x[-1])
plt.ylim(0,1)
plt.legend()
plt.savefig(savedir + 'A_F.png')
plt.close()

x2 = np.linspace(0,4,100)
y2 = F(x2)

plt.plot(x2,y2,color = 'k',label = 'F(x)')
plt.axvline(x = 0.29,color = 'r',label = 'x = 0.29')
plt.ylabel('F(x)')
plt.xlabel('x')
plt.xlim(x2[0],x2[-1])
plt.ylim(0,1)
plt.legend()
plt.savefig(savedir + 'B_F.png')
plt.close()







