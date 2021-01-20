
from ..tool import ch_to_energy,overlap_curve
import numpy as np

class Event(object):

	def __init__(self,ch_n,t,ch):

		self.ch_n = ch_n

		self.dt1 = 2.6e-6   # dead time
		self.dt2 = 10e-6    # the dead time of last channel
		self.t,self.e,self.ch = ch_to_energy(t,ch,self.ch_n['CHANNEL'].values,self.ch_n['E_MIN'].values,self.ch_n['E_MAX'].values)
		self.t_start = self.t.min()
		self.t_stop = self.t.max()

	def __call__(self,binsize = 0.064,bins = None, energy_band = None,channel = False):

		if bins is None:
			bins = np.arange(self.t_start, self.t_stop, binsize)
		size_ = bins[1:] - bins[:-1]

		e_index = np.where(self.ch == self.ch_n['CHANNEL'].values[-1])[0]
		bin_n = np.histogram(self.t[e_index],bins = bins)[0]
		deadtime = bin_n*self.dt2
		e_index = np.where(self.ch != self.ch_n['CHANNEL'].values[-1])[0]
		bin_n = np.histogram(self.t[e_index],bins = bins)[0]
		deadtime = bin_n*self.dt1 + deadtime
		size_ = size_-deadtime
		bin_c = 0.5 * (bins[1:] + bins[:-1])
		if energy_band is not None:
			try:
				if channel:
					e_index = np.where((self.ch>=energy_band[0])&(self.ch<=energy_band[-1]))[0]
				else:
					e_index = np.where((self.e >=energy_band[0])&(self.e<=energy_band[-1]))[0]
				bin_n = np.histogram(self.t[e_index],bins = bins)[0]
				return bin_c,bin_n/size_
			except ValueError:
				rate_list = []
				for e_b in energy_band:
					if channel:
						e_index = np.where((self.ch>=e_b[0])&(self.ch<=e_b[-1]))[0]
					else:
						e_index = np.where((self.e >=e_b[0])&(self.e<=e_b[-1]))[0]
					bin_n = np.histogram(self.t[e_index],bins = bins)[0]
					rate_list.append(bin_n/size_)
				return bin_c,np.vstack(rate_list)
		else:
			bin_n = np.histogram(self.t,bins = bins)[0]
			return bin_c,bin_n/size_

	def ovelap_lightcurve(self,timearray,energy_band = None,channel = False):

		#get die time
		size_ = timearray[:,1]-timearray[:,0]
		e_index = np.where(self.ch == self.ch_n['CHANNEL'].values[-1])[0]
		bin_n = overlap_curve(self.t[e_index],timearray)
		deadtime = bin_n*self.dt2
		e_index = np.where(self.ch != self.ch_n['CHANNEL'].values[-1])[0]
		bin_n = overlap_curve(self.t[e_index],timearray)
		deadtime = bin_n*self.dt1 + deadtime
		size_ = size_-deadtime
		bin_c = 0.5 * (timearray[:,1]+timearray[:,0])
		if energy_band is not None:
			try:
				if channel:
					e_index = np.where((self.ch>=energy_band[0])&(self.ch<=energy_band[-1]))[0]
				else:
					e_index = np.where((self.e >=energy_band[0])&(self.e<=energy_band[-1]))[0]
				bin_n = overlap_curve(self.t[e_index],timearray)
				return bin_c,bin_n/size_
			except ValueError:

				rate_list = []
				for e_b in energy_band:
					if channel:
						e_index = np.where((self.ch>=e_b[0])&(self.ch<=e_b[-1]))[0]
					else:
						e_index = np.where((self.e >=e_b[0])&(self.e<=e_b[-1]))[0]
					bin_n = overlap_curve(self.t[e_index],timearray)
					rate_list.append(bin_n/size_)
				return bin_c,np.vstack(rate_list)
		else:
			bin_n = overlap_curve(self.t,timearray)
			return bin_c,bin_n/size_






