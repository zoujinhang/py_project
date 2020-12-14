
from .clock import Clock
from .background import get_background_f
from .satellite import Detectors,Geometry,Locate
#from .perception import Event
from .perception import Event,try_to_trig,bayesian_trig
import numpy as np



class Search(object):


	def __init__(self,data,pd_position_data):

		#self.data = data
		self.name = np.array(list(data.keys()))
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
		self.energy_band = [[5,900],[5,50],[44,300]] #[300,900]
		self.candidate = []
		self.time_list = []

		print('search 1')
		for index ,(start,stop) in enumerate(self.window_list):
			print('search 1 ---',index)
			bins_baseline = np.arange(start,stop,self.binsize_baseline)
			bins_baseline_c = 0.5 * (bins_baseline[1:] + bins_baseline[:-1])
			bb_cc = np.concatenate(([start], bins_baseline_c, [stop]))
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
					rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
					bs_f = get_background_f(bb_cc,rate_b)
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

			#SNR_arr = np.vstack(SNR_arr)
			#lc_arr = np.vstack(lc_arr)
			#bs_arr = np.vstack(bs_arr)
			self.candidate.append([SNR_arr,lc_arr,bs_arr])
			self.time_list.append(bins_lightcurve_c)
		print('search 1 end!')

	def locate(self):

		#rate = np.zeros((len(self.energy_band),len(self.name)))
		#bs = np.zeros((len(self.energy_band),len(self.name)))

		location5_50 = []
		location50_300 = []

		for index,lc_t in enumerate(self.time_list):
			print('locat ',index)
			edges_ = self.edges_list[index]
			#during = edges_[-1] - edges_[0]
			index_ = np.where((lc_t>=edges_[0]-1)&(lc_t<=edges_[-1]+1))[0]
			SNR_arr,lc_arr,bs_arr = self.candidate[index]
			#e = SNR_arr[0]
			index_list5_50 = []
			index_snr5_50 = None
			index_list50_300 = []
			index_snr50_300 = None
			for i in index_:
				#print('5-50',SNR_arr[1][i])
				index_snr = np.where(SNR_arr[1][i]>=4.)[0]
				if len(index_snr)>=3 and len(index_list5_50)*self.binsize_lightcurve <= 2:
					index_list5_50.append(i)
					if index_snr5_50 is None:
						index_snr5_50 = index_snr
					else:
						if len(index_snr)>len(index_snr5_50):
							index_snr5_50=index_snr
				#print('50-300',SNR_arr[2][i])
				index_snr = np.where(SNR_arr[2][i]>=4.)[0]
				if len(index_snr)>=3 and len(index_list50_300)*self.binsize_lightcurve <= 2: #good for location of energy band 50-300
					index_list50_300.append(i)
					if index_snr50_300 is None:
						index_snr50_300 = index_snr
					else:
						if len(index_snr)>len(index_snr50_300):
							index_snr50_300=index_snr
			if len(index_list5_50)>0:

				lc_5_50 = lc_arr[1][index_list5_50]
				m_rate_5_50 = lc_5_50.mean(axis = 0)
				bs_5_50 = bs_arr[1][index_list5_50]
				m_bs_5_50 = bs_5_50.mean(axis = 0)
				t = (lc_t[index_list5_50]).mean()
				during = len(index_list5_50)*self.binsize_lightcurve
				detectorlist = self.name[index_snr5_50]
				location = self.location.locate(t,m_rate_5_50,m_bs_5_50,0,self.energy_band[1],during=during,detector_list=detectorlist)
				location5_50.append(location)
			else:
				location5_50.append(None)

			if len(index_list50_300)>0:

				lc_50_300 = lc_arr[2][index_list50_300]
				m_rate_50_300 = lc_50_300.mean(axis = 0)
				bs_50_300 = bs_arr[2][index_list50_300]
				m_bs_50_300 = bs_50_300.mean(axis = 0)
				t = (lc_t[index_list50_300]).mean()
				during = len(index_list5_50)*self.binsize_lightcurve
				detectorlist = self.name[index_snr50_300]
				location = self.location.locate(t,m_rate_50_300,m_bs_50_300,1,self.energy_band[2],during=during,detector_list=detectorlist)
				location50_300.append(location)
			else:
				location50_300.append(None)

		return location5_50,location50_300



