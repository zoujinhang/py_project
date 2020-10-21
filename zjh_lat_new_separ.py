
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import os
import Data_analysis.file as myfile
matplotlib.use('Agg')
from multiprocessing import Pool
import sys
'''
def run(a):
	print(a)
	plt.plot(a,0)
	plt.savefig('/home/laojin/my_lat/'+str(a)+'.png')
	plt.close()

a = [1,2,3,5]
pool = Pool(2)
pool.map(run,a)
pool.close()
pool.join()
'''
paths = sys.path
c_lib_link = None
for path in paths:
	fand = path+'/Data_analysis/Separate_source/WT_c/'
	print(fand)
	if os.path.exists(fand):
		sonamelist = myfile.findfile(fand,'WT_weight.so*')
		print(sonamelist)
		if len(sonamelist)>0:

			c_lib_link = fand+sonamelist[0]
			print('the C lib link is ',c_lib_link)
			break
if c_lib_link is not None:
	print(c_lib_link)
else:
	print('can not find the C lib of WT_weight.os!')







