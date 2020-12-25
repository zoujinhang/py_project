
from ..tool import ch_to_energy
import numpy as np

class Event(object):

	def __init__(self,ch_n,t,ch):

		self.ch_n = ch_n

		self.dt1 = 2.6e-6  # dead time
		self.dt2 = 10e-6   #the dead time of last channel
		self.t,self.e,self.ch = ch_to_energy(t,ch,self.ch_n['CHANNEL'].values,self.ch_n['E_MIN'].values,self.ch_n['E_MAX'].values)
		self.t_start = self.t.min()
		self.t_stop = self.t.max()

	def __call__(self,binsize = 0.064,bins = None, energy_band = None,channel = False):

		if bins is None:
			bins = np.arange(self.t_start, self.t_stop, binsize)
		size_ = bins[1:] - bins[:-1]

		e_index = np.where(self.ch == self.ch_n['CHANNEL'].values[-1])[0]
		t = self.t[e_index]
		bin_n = np.histogram(t,bins = bins)[0]
		deadtime1 = bin_n*self.dt2
		e_index = np.where(self.ch != self.ch_n['CHANNEL'].values[-1])[0]
		t = self.t[e_index]
		bin_n = np.histogram(t,bins = bins)[0]
		deadtime = bin_n*self.dt1 + deadtime1
		size_ = size_-deadtime

		bin_c = 0.5 * (bins[1:] + bins[:-1])
		t = self.t
		if energy_band is not None:
			if channel:
				e_index = np.where((self.ch>=energy_band[0])&(self.ch<=energy_band[-1]))[0]
			else:
				e_index = np.where((self.e >=energy_band[0])&(self.e<=energy_band[-1]))[0]
			t = t[e_index]

		bin_n = np.histogram(t,bins = bins)[0]
		rate = bin_n/size_
		return bin_c,rate


















