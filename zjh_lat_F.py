import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import gamma

def K(s):
	return mp.besselk(5/3,s)#贝塞尔函数

def F(x):
	f = np.zeros(len(x))
	for index,value in enumerate(x):
		f[index] = value*mp.quad(K,[value,np.inf])#求积分
	return f

def F_lt_1(x):
	return 4*np.pi*(x/2)**(1/3)/np.sqrt(3)/gamma(1/3)#远小于1近似公式
def F_gt_1(x):
	return (np.pi/2)**(1/2)*np.exp(-x)*x**(1/2)#远大于1近似公式

savedir = '/home/laojin/my_lat/synchrotron_radiation/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

#a  = mp.besselk(5/3,0)
#print(a)

x = np.logspace(-3,1,100)
y = F(x)

plt.plot(x,y,color = 'k',label = 'F(x)')
plt.axvline(x = 0.3,color = 'r',label = 'x = 0.3')
plt.plot(x[x>1],F_gt_1(x[x>1]),color = 'g',label = 'Approximation of F(x) with x>>1')
plt.plot(x[x<1],F_lt_1(x[x<1]),color = 'y',label = 'Approximation of F(x) with x<<1')
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







