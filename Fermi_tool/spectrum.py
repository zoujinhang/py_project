from astropy.io import fits
import numpy as np

class Rsp(object):
	'''
	
	'''
	def __init__(self,rsp_link):
		hl = fits.open(rsp_link)
		matr = hl[2].data['MATRIX']
		matr[-1] = np.zeros(128)
		matr = np.vstack(matr)
		self.matr = np.mat(matr)
		self.e_lo = hl[2].data['ENERG_LO']
		self.e_hi = hl[2].data['ENERG_HI']
		self.e_min = hl[1].data['E_MIN']
		self.e_max = hl[1].data['E_MAX']
	
	def transform(self,model):
		model_mat = np.mat(model)
		T_model = model_mat*self.matr
		return np.array(T_model)[0]

	









