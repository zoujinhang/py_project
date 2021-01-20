

import numpy as np
from scipy.sparse import csc_matrix,eye,diags
from scipy.sparse.linalg import spsolve
import sys
import ctypes
import os
from ..file import findfile

paths = sys.path
c_lib_link = None
for path in paths:
	fand = path+'/daily_search/background/c_lib/'
	if os.path.exists(fand):
		sonamelist = findfile(fand,'c_baseline.so.6')
		if len(sonamelist)>0:

			c_lib_link = fand+sonamelist[0]
			print('the C lib link is ',c_lib_link)
			break
if c_lib_link is not None:
	clib = ctypes.cdll.LoadLibrary(c_lib_link)
else:
	print('can not find the C lib of c_baseline.so.6!')


def TD_baseline(time,rate,lam = None,hwi = None,it = None,inti = None,lambda_2 = None):
	'''

	:param time:
	:param rate:
	:param lam:
	:param hwi:
	:param it:
	:param inti:
	:return:
	'''
	dt = time[1]-time[0]
	during = np.max(time)-np.min(time)
	if(lam is None):
		lam = 0.1/dt**1.5
	else:
		lam = lam/dt**1.5

	if(lambda_2 is None):
		lambda_2 = 2/dt**1.5
	else:
		lambda_2 = lambda_2/dt**1.5

	if(it is None):
		it = 14
	if(inti is None):

		fillpeak_int = int(during/2)
		inti = 2
	else:
		fillpeak_int = int(during/inti)

	if (hwi is None):
		hwi = int(250 / inti)
	else:
		hwi = int(hwi / inti)

	if fillpeak_int > len(rate):
		fillpeak_int = len(rate)
	# if len(rate)<50 and len(rate)>5:
	#	fillpeak_int = 5
	if len(rate) <= 5:
		fillpeak_int = len(rate)

	if(lam < 1):
		lam = 1
	if lambda_2 <1:
		lambda_2 = 1
	return baseline_kernel(rate,lambda_=lam,hwi=hwi,it = it,int_ = fillpeak_int,lambda_2=lambda_2)


def get_smooth(spectra,lambda_):
	'''

	:param spectra:
	:param lambda_:
	:return:
	'''
	spectra = np.array(spectra)
	m = spectra.shape[0]
	w = np.ones(m)
	smooth = WhittakerSmooth(spectra,w,lambda_)
	cs = spectra-smooth
	cs_mean = cs.mean()
	cs_std = cs.std()
	for i in range(3):
		cs_index = np.where((cs>cs_mean+(1+1*i)*cs_std)|(cs<cs_mean-(1+1*i)*cs_std))
		w[cs_index] = 0
		smooth = WhittakerSmooth(spectra,w,lambda_)
		cs = spectra-smooth
		cs_mean = cs[w!=0].mean()
		cs_std = cs[w!=0].std()
	return smooth


def baseline_kernel(spectra,lambda_,hwi,it,int_,lambda_2):
	'''

	:param spectra:
	:param lambda_:
	:param hwi:
	:param it:
	:param int_:
	:return:
	'''
	spectra = np.array(spectra)
	w = np.ones(spectra.shape[0])
	spectra = WhittakerSmooth(spectra,w,lambda_)

	if it != 1 :
		d1 = np.log10(hwi)
		d2 = 0
		w = np.ceil(np.concatenate((10**(d1+np.arange(0,it-1,1)*(d2-d1)/(np.floor(it)-1)),[d2])))
		w = np.array(w,dtype = int)
	else:
		w = np.array([hwi],dtype = int)
	#print(w)

	lims = np.linspace(0,spectra.size -1,int_+1)
	lefts = np.array(np.ceil(lims[:-1]),dtype = int)#This is the index value
	rights = np.array(np.floor(lims[1:]),dtype = int)#Same as above
	minip = (lefts+rights)*0.5#The index

	if clib is not None:
		c_spectra = (ctypes.c_double * spectra.size)(*list(spectra))
		len_spec = (ctypes.c_int32)(spectra.size)
		c_lefts = (ctypes.c_int32 * int_)(*list(lefts))
		c_rights = (ctypes.c_int32 * int_)(*list(rights))
		c_xx = (ctypes.c_double * int_)()
		c_w =  (ctypes.c_int32 * it)(*list(w))
		len_w = (ctypes.c_int32)(it)
		len_xx = (ctypes.c_int32)(int_)
		clib.c_baseline_kernel(c_spectra,len_spec,c_xx,c_lefts,c_rights,len_xx,c_w,len_w)
		xx = np.array(c_xx)
	else:

		xx = np.zeros(int_)
		for i in range(int_):
			xx[i] = spectra[lefts[i]:rights[i] + 1].mean()
		for i in range(it):
			# Current window width
			w0 = w[i]
			# Point-wise iteration to the right
			for j in range(1,int_-1):
				# Interval cut-off close to edges
				v = min([j,w0,int_-j-1])
				# Baseline suppression
				a = xx[j-v:j+v+1].mean()
				xx[j] = min([a,xx[j]])
			for j in range(1,int_-1):
				k = int_-j-1
				v = min([j,w0,int_-j-1])
				a = xx[k-v:k+v+1].mean()
				xx[k] = min([a,xx[k]])

	minip = np.concatenate(([0],minip,[spectra.size-1]))
	xx = np.concatenate((xx[:1],xx,xx[-1:]))
	index = np.arange(0,spectra.size,1)
	xxx = np.interp(index,minip,xx)
	w = np.ones(xxx.shape[0])
	return WhittakerSmooth(xxx, w, lambda_2)

def WhittakerSmooth(x,w,lambda_):
	'''

	:param x: array
	:param w: array .An array of weights corresponding to the values
	:param lambda_: Smoothing parameter
	:return: array Smoothing results
	'''

	X=np.mat(x)
	m=X.size
	#i=np.arange(0,m)
	E=eye(m,format='csc')
	D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
	W=diags(w,0,shape=(m,m))
	A=csc_matrix(W+(lambda_*D.T*D))
	B=csc_matrix(W*X.T)
	background=spsolve(A,B)

	return np.array(background)







