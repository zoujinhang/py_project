

import Data_analysis.file as myfile
import sys
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt

paths = sys.path
c_lib_link = None
for path in paths:
	fand = path+'/Data_analysis/Separate_source/WT_c/'
	if os.path.exists(fand):
		sonamelist = myfile.findfile(fand,'WT_weight.so*')
		if len(sonamelist)>0:

			c_lib_link = fand+sonamelist[0]
			print('the C lib link is ',c_lib_link)
			break
if c_lib_link is not None:
	clib = ctypes.cdll.LoadLibrary(c_lib_link)
else:
	print('can not find the C lib of WT_weight.os!')


t = np.arange(0,10,0.5)
dt_ = t[1:]-t[:-1]

print('t',t)
print('dt',dt_)

len_t = len(t)
len_dt = len(dt_)
t = (ctypes.c_double * len_t)(*list(t))
dt = (ctypes.c_double * len_dt)(*list(dt_))
wt = (ctypes.c_double * len_t)()
sigma = 0.5/(1)**2
#clib.work_WT(wt,(ctypes.c_double)(2.0),dt,t,len_t,len_dt)
clib.work_WT(wt,(ctypes.c_double)(sigma),dt,t,len_t,len_dt)
print('wt',*list(wt))

plt.plot(t,list(wt))
plt.savefig('wt_check.png')
plt.close()

