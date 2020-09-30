
import numpy as np
from scipy.sparse import csc_matrix,eye,diags
from scipy.sparse.linalg import spsolve
import operator
from scipy import stats
from astropy.stats import sigma_clip,mad_std
import os
import re


def transpose(data):
	da = [[row[col] for row in data] for col in range(len(data[0]))]
	return da

def readcol(file_name):
	'''
	read txt
	:param file_name:
	:return:
	'''
	c = []
	f = open(file_name,'r')
	for line in f:
		a = line
		b = [i for i in re.split('[\t,\n,\s]',a) if i != '']
		#print(b)
		c = c + [b]
	f.close()
	c = [i for i in c if i != []]
	data_all = transpose(c)
	nl = len(data_all)
	for k in range(nl):
		for l in range(len(data_all[k])):
			try:
				data_all[k][l] = int(data_all[k][l])
			except ValueError:
				try:
					data_all[k][l] = float(data_all[k][l])
				except ValueError:
					data_all[k][l] = data_all[k][l]
	return data_all


def get_energy_of_ch(time,e1,e2):
	'''
	
	:param time: time
	:param e1:
	:param e2:
	:return:
	'''
	numb = len(time)
	energy_random_arr = np.random.random_sample(numb)
	energy_array = e1 + (e2-e1)*energy_random_arr
	return energy_array

def ch_to_energy(time,ch,ch_n,e1,e2):
	'''
	
	:param time:
	:param ch:
	:param ch_n:
	:param e1:
	:param e2:
	:return:
	'''
	new_t = np.array([])
	new_energy = np.array([])
	for index,channel in enumerate(ch_n):
		ch_t_index = np.where(ch == channel)
		ch_t = time[ch_t_index]
		energy_array = get_energy_of_ch(ch_t,e1[index],e2[index])
		new_t = np.concatenate((new_t,ch_t))
		new_energy = np.concatenate((new_energy,energy_array))
	index_all = np.argsort(new_t)
	new_t = new_t[index_all]
	new_energy = new_energy[index_all]
	return new_t,new_energy


def findfile(dir1, feature):
	'''

	:param dir1:
	:param feature:
	:return:
	'''
	if (os.path.exists(dir1)):
		dirnames = os.listdir(dir1)
		filelist = []
		fil_number = 0
		fil_result_number = 0
		featurelist = [i for i in re.split('[*]', feature) if i != '']
		for_number = len(featurelist)
		fileresult = [[] for i in range(for_number)]
		for eve in range(for_number):
			if (eve == 0):
				fe_number = len(featurelist[eve])
				for sample in dirnames:
					if (os.path.isfile(dir1 + sample)):
						filelist.append(sample)
						fil_number = fil_number + 1
				if (fil_number != 0):
					for i in filelist:
						i_number = len(i)
						n = i_number - fe_number + 1
						for j in range(n):
							if (i[j:j + fe_number] == featurelist[eve]):
								fileresult[eve].append(i)
								fil_result_number = fil_result_number + 1
								break
					# print('1----------',fileresult[eve])#------------------------
					if (fil_result_number == 0):
						print(
							'we do not find any file that has the feature with [' + feature + ']!\n')
						return []
					else:
						fil_result_number = 0
				else:
					print('there is no file in this dir ! \n')
					return []
			else:
				fe_number = len(featurelist[eve])
				for i in fileresult[eve - 1]:
					i_number = len(i)
					n = i_number - fe_number + 1
					for j in range(n):
						if (i[j:j + fe_number] == featurelist[eve]):
							fileresult[eve].append(i)
							fil_result_number = fil_result_number + 1
							break
				if (fil_result_number == 0):
					print('we do not find any file that has the feature with [' + feature + ']!\n')
					return []
				else:
					fil_result_number = 0
		return fileresult[for_number - 1]
	else:
		print('do not find the dir named [' + dir1 + ']!\n')
		return False


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

def get_w(cs,sigma):
	return np.exp(-0.5/(sigma**2)*(cs)**2)

def TD_bs(t,rate,it_ = 1,lambda_=4000,sigma = False,hwi = None,it = None,inti = None):
	dt = t[1]-t[0]
	t_c,cs,bs = TD_baseline(t,rate,hwi = hwi,it = it ,inti =inti)
	mask = sigma_clip(cs, sigma=5, maxiters=5, stdfunc=mad_std).mask
	myfilter = list(map(operator.not_, mask))
	lc_median_part = cs[myfilter]
	loc, scale = stats.norm.fit(lc_median_part)
	for i in range(it_):
		w = get_w(cs,scale)
		bs = WhittakerSmooth(rate,w,lambda_=lambda_/dt**1.5)
		cs = rate - bs
		mask = sigma_clip(cs, sigma=5, maxiters=5, stdfunc=mad_std).mask
		myfilter = list(map(operator.not_, mask))
		lc_median_part = cs[myfilter]
		loc, scale = stats.norm.fit(lc_median_part)
		
	if sigma:
		return cs,bs,scale
	else:
		return cs,bs


def TD_baseline(time,rate,lam = None,hwi = None,it = None,inti = None):
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
	if(lam is None):
		lam = 100/dt**1.5
	else:
		lam = lam/dt**1.5
	if(hwi is None):
		hwi = int(20/dt)
	else:
		hwi = int(hwi/dt)
	if(it is None):
		it = 5
	if(inti is None):

		fillpeak_int = int(len(rate)/10)

	else:
		fillpeak_int =inti
	if(lam < 1):
		lam = 1
	bs = baseline_kernel(rate,lambda_=lam,hwi=hwi,it = it,int_ = fillpeak_int)
	return time,rate-bs,bs


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


def baseline_kernel(spectra,lambda_,hwi,it,int_):
	'''
	
	:param spectra:
	:param lambda_:
	:param hwi:
	:param it:
	:param int_:
	:return:
	'''
	spectra = np.array(spectra)
	spectra = get_smooth(spectra,lambda_)

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
	xx = np.zeros(int_)
	for i in range(int_):
		xx[i] = spectra[lefts[i]:rights[i]+1].mean()

	
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
	return xxx



def get_band(band ,numb ,ovelap=0.5 ,scale='log'):
	
	if scale == 'log':
		return get_band_in_log(band ,numb ,ovelap)
	else:
		return get_band_in_unif(band ,numb ,ovelap)


def get_band_in_log(band ,numb ,ovelap=0.5):
	if ovelap >= 0.95:
		ovelap = 0.95
	non_ovelap = 1- ovelap
	log_band = np.log10(band)
	len_ = np.abs(log_band[-1] - log_band[0])
	x = len_ / (non_ovelap * (numb - 1) + 1)
	log_el_edges = np.linspace(log_band[0], log_band[-1] - x, numb)
	log_eh_edges = log_el_edges + x
	return np.vstack([10 ** log_el_edges, 10 ** log_eh_edges]).T


def get_band_in_unif(band, numb, ovelap=0.5):
	if ovelap >= 0.95:
		ovelap = 0.95
	non_ovelap = 1 - ovelap
	len_ = np.abs(band[-1] - band[0])
	x = len_ / (non_ovelap * (numb - 1) + 1)
	el_edges = np.linspace(band[0], band[-1] - x, numb)
	eh_edges = el_edges + x
	
	return np.vstack([el_edges, eh_edges]).T





