
from .clock import Clock
from .background import get_background_f
from .satellite import Detectors,Geometry,Locate
from .perception import Event
from .perception import Event,try_to_trig,bayesian_trig
import numpy as np



class Search(object):


	def __init__(self,data,pd_position_data):

		#self.data = data
		self.name = data.keys()
		self.geometry = Geometry(pd_position_data)
		self.location = Locate(self.geometry)
		self.detector = self.geometry.detectors

		self.lc_window_list = try_to_trig(data,self.detector)
		print('try_to_trig :\n',self.lc_window_list)
		self.edges_list,self.window_list,self.name_list,self.lc_wind_index_list = bayesian_trig(data,self.lc_window_list,self.detector)

		self.sigma = 4
		self.strong = 7
		self.pfstrong = 10
		self.wt = 0.1
		self.binsize_baseline = 1
		self.binsize_lightcurve = 0.05
		self.energy_band = [[5,900],[5,50],[50,300],[300,900]]
		self.candidate = []
		self.time_list = []
		for index ,(start,stop) in enumerate(self.window_list):

			bins_baseline = np.arange(start,stop,self.binsize_baseline)
			bins_lightcurve = np.arange(start,stop,self.binsize_lightcurve)
			bins_lightcurve_c = 0.5*(bins_lightcurve[1:]+bins_lightcurve[:-1])
			SNR_list = []
			lc_list = []
			bs_list = []

			for detei in self.name:

				ni = data[detei]
				ch_E = ni['ch_E']
				t = ni['events']['TIME'].values
				ch = ni['events']['PHA'].values
				lightcurve = Event(ch_E,t,ch)
				detector_SNR_list = []
				detector_lc_list = []
				detector_bs_list = []
				for e_band in self.energy_band:
					rate_b = lightcurve(bins = bins_baseline,energy_band=e_band)[1]
					rate_b = np.concatenate((rate_b[:1],rate_b))
					bs_f = get_background_f(bins_baseline,rate_b)
					t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=e_band)
					bs = bs_f(t_cc)
					scale = np.sqrt(bs/self.binsize_lightcurve)
					detector_SNR_list.append((rate_c-bs)/scale)
					detector_bs_list.append(bs)
					detector_lc_list.append(rate_c)

				SNR_list.append(detector_SNR_list)
				lc_list.append(detector_lc_list)
				bs_list.append(detector_bs_list)

			SNR_arr = []
			lc_arr = []
			bs_arr = []
			for i in range(len(self.energy_band)):

				e_band_SNR = []
				e_band_bs = []
				e_band_lc = []

				for j in range(len(self.name)):

					e_band_SNR.append(SNR_list[j][i])
					e_band_bs.append(bs_list[j][i])
					e_band_lc.append(lc_list[j][i])

				e_band_SNR = np.vstack(e_band_SNR).T
				e_band_bs = np.vstack(e_band_bs).T
				e_band_lc = np.stack(e_band_lc).T

				SNR_arr.append(e_band_SNR)
				lc_arr.append(e_band_lc)
				bs_arr.append(e_band_bs)

			SNR_arr = np.vstack(SNR_arr)
			lc_arr = np.vstack(lc_arr)
			bs_arr = np.vstack(bs_arr)
			self.candidate.append([SNR_arr,lc_arr,bs_arr])
			self.time_list.append(bins_lightcurve_c)





					



