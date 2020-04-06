'''

'''
import numpy as np
from ..Baseline import WhittakerSmooth

class WT_kernel(object):

	def __init__(self,t,dt = 1.024):

		self.dt = dt
		self.t = np.array(t)
		self.wt = self.work_WT()
		rate = 1.0/self.wt
		w = np.ones(rate.size)
		self.rate = WhittakerSmooth(rate,w,2)

	def get_wt(self):
		return self.wt

	def get_rate(self):
		return self.rate

	def weight(self,x,x0,sigma):
		return np.exp(-(x-x0)**2/(2*sigma**2))

	def work_WT(self):
		#这里的卷积有些慢。
		dt = self.t[1:]-self.t[:-1]
		wt = []
		for index,value in enumerate(self.t):
			ww = self.weight(self.t,value,self.dt)  #长度为t的长度。
			dti = np.insert(dt,index,0)
			ww_sum = ww.sum()-1
			wti = (dti*ww).sum()/ww_sum
			wt.append(wti)
		wt = np.array(wt)
		return wt