
#from .clock import Clock
from .background import TD_baseline
from .satellite import Geometry,Locate,Sky_map
#from .perception import Event
from .perception import Event,try_to_trig,bayesian_trig
import numpy as np
from astropy.coordinates import SkyCoord
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



class Search(object):


	def __init__(self,data,pd_position_data,marker = None):

		#self.data = data
		self.marker = marker
		self.name = np.array(list(data.keys()))
		self.geometry = Geometry(pd_position_data)
		self.location = Locate(self.geometry)
		self.detector = self.geometry.detectors
		self.clock = self.geometry.clock
		self.lc_window_list = try_to_trig(data,self.detector)
		if self.marker is not None:
			print(self.marker + ' try_to_trig number:',len(self.lc_window_list))
		self.edges_list,self.window_list,self.name_list,self.lc_wind_index_list = bayesian_trig(data,self.lc_window_list,self.detector)
		if self.marker is not None:
			print(self.marker + ' bayesian trig number:',len(self.edges_list))
		self.sigma = 4
		self.strong = 7
		self.pfstrong = 10
		self.wt = 0.1
		#self.binsize_baseline = 1
		self.binsize_lightcurve = 0.05
		self.energy_band = [[8,105],[8,32],[33,84]] #[300,900]
		self.cenergies = [[5,50],[44,300]]
		self.candidate = []
		self.time_list = []

		#print('search 1')
		for index ,(start,stop) in enumerate(self.window_list):
			if self.marker is not None:
				print(self.marker + ' search ---',index)
			#bins_baseline = np.arange(start,stop,self.binsize_baseline)
			#bins_baseline_c = 0.5 * (bins_baseline[1:] + bins_baseline[:-1])
			#bb_cc = np.concatenate(([start], bins_baseline_c, [stop]))
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
					#rate_b = lightcurve(bins = bins_baseline,energy_band=e_band,channel = True)[1]
					#rate_b = np.concatenate((rate_b[:1],rate_b,rate_b[-1:]))
					#bs_f = get_background_f(bb_cc,rate_b)
					t_cc,rate_c = lightcurve(bins = bins_lightcurve,energy_band=e_band,channel = True)
					#bs = bs_f(t_cc)
					bs = TD_baseline(t_cc,rate_c)
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
		#print('search 1 end!')

		self.location5_50,self.location50_300 = self.locate()

	def track_one(self,trg,name = None):

		return_lc_t = []
		location5_50 = []
		location50_300 = []
		edges_list = []
		new_namelist = []
		candidate = []
		for index,lc_t in enumerate(self.time_list):
			edges_ = self.edges_list[index]
			namelist = self.name_list[index]
			location = False
			if location == False:
				if self.location5_50[index] is not None:
					ra_,dec_,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
					po_rcat = SkyCoord(ra = ra_rcat,dec = dec_rcat,frame = 'icrs',unit = 'deg')
					sep = po_rcat.separation(trg).deg
					index1 = np.argmin(sep)
					d_chi2 = chi2 - chi2.min()
					if d_chi2[index1] <= 9.21 and sep[index1] <=2:
						location = True
			elif location == False:
				if self.location50_300[index] is not None:
					ra_,dec_,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
					po_rcat = SkyCoord(ra = ra_rcat,dec = dec_rcat,frame = 'icrs',unit = 'deg')
					sep = po_rcat.separation(trg).deg
					index1 = np.argmin(sep)
					d_chi2 = chi2 - chi2.min()
					if d_chi2[index1] <= 9.21 and sep[index1] <=2:
						location = True
			elif location == False:
				seq = self.geometry.get_separation_with_time(edges_[0], trg)
				earth_p = self.geometry.get_earth_point(edges_[0])
				n_ni = len(namelist)
				earth_ra,earth_dec,radius_deg = earth_p[1:]
				xyz_position = SkyCoord(ra = earth_ra,dec = earth_dec,frame = 'icrs',unit = 'deg')
				seq_with_earth = xyz_position.separation(trg).deg
				if seq_with_earth>radius_deg:
					seqi = seq.iloc[0]
					seqi_value = seqi.values[1:]
					seqi_value_min = np.min(seqi_value)+4
					seq_values = seqi[namelist].values
					ne = seq_values[seq_values<65]
					if (len(ne)>=0.333*n_ni):
						if (np.min(ne)<=seqi_value_min):
							location = True
			if location:
				return_lc_t.append(lc_t)
				location5_50.append(self.location5_50[index])
				location50_300.append(self.location50_300[index])
				edges_list.append(edges_)
				new_namelist.append(namelist)
				candidate.append(self.candidate[index])
		return Target(self.name,return_lc_t,candidate,edges_list,new_namelist,location5_50,location50_300,target = trg,name = name,geometry=self.geometry)


	def locate(self):

		#rate = np.zeros((len(self.energy_band),len(self.name)))
		#bs = np.zeros((len(self.energy_band),len(self.name)))

		location5_50 = []
		location50_300 = []

		for index,lc_t in enumerate(self.time_list):
			#print('locat ',index)
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
				index_snr = np.where(SNR_arr[1][i]>=4.5)[0]
				if len(index_snr)>=3 and len(index_list5_50)*self.binsize_lightcurve <= 2:
					index_list5_50.append(i)
					if index_snr5_50 is None:
						index_snr5_50 = index_snr
					else:
						if len(index_snr)>len(index_snr5_50):
							index_snr5_50=index_snr
				#print('50-300',SNR_arr[2][i])
				index_snr = np.where(SNR_arr[2][i]>=4.5)[0]
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
				location = self.location.locate(t,m_rate_5_50,m_bs_5_50,0,self.cenergies[0],detector_list=detectorlist,during=during)
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
				location = self.location.locate(t,m_rate_50_300,m_bs_50_300,1,self.cenergies[1],detector_list=detectorlist,during=during)
				location50_300.append(location)
			else:
				location50_300.append(None)

		return location5_50,location50_300

	def save_triggers(self,savename):

		trig_time = []
		for i in self.edges_list:
			trig_time.append(i[0])
		trig_time = np.array(trig_time)
		utc_time = self.clock.met_to_utc(trig_time).fits
		c = {'met':trig_time,'utc':utc_time,'trig_detector':self.name_list}
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename,index=False)

	def save_location_5_50_chi2(self,index,savename):

		if self.location5_50[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
			c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def save_location_5_50(self,index,savename):

		if self.location5_50[index] is not None:

			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
			c = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)


	def save_location_50_300_chi2(self,index,savename):

		if self.location50_300[index] is not None:

			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
			c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)


	def save_location_50_300(self,index,savename):

		if self.location50_300[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
			c = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def get_candidate_number(self):

		return len(self.time_list)

	def get_location_5_50(self,index):
		return self.location5_50[index]

	def get_location_50_300(self,index):

		return self.location50_300[index]

	def get_lc_time(self,index):

		return self.time_list[index]

	def get_lc_rate_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_5_50 = lc_arr[1]
		return lc_5_50.T

	def get_lc_bs_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_5_50 = bs_arr[1]
		return bs_5_50.T

	def get_SNR_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_5_50 = SNR_arr[1]
		return SNR_5_50.T

	def save_candidate_5_50_signal(self,index,savename):

		c = {'met':self.time_list[index],'utc':self.clock.met_to_utc(self.time_list[index]).fits}
		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_ = SNR_arr[1]
		bs_ = bs_arr[1]
		lc_ = lc_arr[1]
		for index,ni in enumerate(self.name):
			c[ni+'_lc'] = lc_[:,index]
			c[ni+'_bs'] = bs_[:,index]
			c[ni+'_SNR'] = SNR_[:,index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename,index=False)

	def get_lc_rate_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_50_300 = lc_arr[2]
		return lc_50_300.T

	def get_lc_bs_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_50_300 = bs_arr[2]
		return bs_50_300.T

	def get_SNR_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_50_300 = SNR_arr[2]
		return SNR_50_300.T

	def save_candidate_50_300_signal(self, index, savename):

		c = {'met': self.time_list[index], 'utc': self.clock.met_to_utc(self.time_list[index]).fits}
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		SNR_ = SNR_arr[2]
		bs_ = bs_arr[2]
		lc_ = lc_arr[2]
		for index, ni in enumerate(self.name):
			c[ni + '_lc'] = lc_[:, index]
			c[ni + '_bs'] = bs_[:, index]
			c[ni + '_SNR'] = SNR_[:, index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename, index=False)


	def get_lc_rate_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_ = lc_arr[0]
		return lc_.T

	def get_lc_bs_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_ = bs_arr[0]
		return bs_.T

	def get_SNR_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_ = SNR_arr[0]
		return SNR_.T

	def save_candidate_signal(self, index, savename):

		c = {'met': self.time_list[index], 'utc': self.clock.met_to_utc(self.time_list[index]).fits}
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		SNR_ = SNR_arr[0]
		bs_ = bs_arr[0]
		lc_ = lc_arr[0]
		for index, ni in enumerate(self.name):
			c[ni + '_lc'] = lc_[:, index]
			c[ni + '_bs'] = bs_[:, index]
			c[ni + '_SNR'] = SNR_[:, index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename, index=False)

	def plot_candidate_50_300_lc(self,index,savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_arri = lc_arr[2].T
		bs_arri = bs_arr[2].T
		#SNR_arri = SNR_arr[2].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(lc_arri)):
			detec = self.name[j]
			plt.subplot(6,2,j+1)
			plt.plot(lc_t-lc_t[0],lc_arri[j],label = 'lc '+detec)
			plt.plot(lc_t-lc_t[0],bs_arri[j],label = 'bs '+detec)
			if detec in namelist:
				plt.axvline(x = edges_[0]-lc_t[0],color = 'r')
				plt.axvline(x = edges_[-1]-lc_t[0],color = 'g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_50_300_SNR(self,index,savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		#lc_arri = lc_arr[2].T
		#bs_arri = bs_arr[2].T
		SNR_arri = SNR_arr[2].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.name[j]
			plt.subplot(6,2,j+1)
			plt.plot(lc_t-lc_t[0],SNR_arri[j],label = 'SNR '+detec)
			if detec in namelist:
				plt.axvline(x = edges_[0]-lc_t[0],color = 'r')
				plt.axvline(x = edges_[-1]-lc_t[0],color = 'g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_5_50_lc(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		lc_arri = lc_arr[1].T
		bs_arri = bs_arr[1].T
		# SNR_arri = SNR_arr[1].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(lc_arri)):
			detec = self.name[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], lc_arri[j], label='lc ' + detec)
			plt.plot(lc_t - lc_t[0], bs_arri[j], label='bs ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_5_50_SNR(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		# lc_arri = lc_arr[1].T
		# bs_arri = bs_arr[1].T
		SNR_arri = SNR_arr[1].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.name[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], SNR_arri[j], label='SNR ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_lc(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		lc_arri = lc_arr[0].T
		bs_arri = bs_arr[0].T
		# SNR_arri = SNR_arr[2].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(lc_arri)):
			detec = self.name[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], lc_arri[j], label='lc ' + detec)
			plt.plot(lc_t - lc_t[0], bs_arri[j], label='bs ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_SNR(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		# lc_arri = lc_arr[0].T
		# bs_arri = bs_arr[0].T
		SNR_arri = SNR_arr[0].T
		lc_t = self.time_list[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig ' + utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.name[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], SNR_arri[j], label='SNR ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_sky_map(self,index,savename):

		edges_ = self.edges_list[index]
		namelist = self.name_list[index]
		smp = Sky_map(figsize=(10,10))
		smp.add_subplot(2,1,1)
		utc = self.clock.met_to_utc(edges_[0])
		smp.title('Energy band 5-50 kev, trigger time ' + utc.fits)
		if self.location5_50[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
			d_chi2 = chi2-chi2.min()
			smp.tricontour(ra_rcat,dec_rcat,d_chi2,[0,2.3,4.61,9.21])
			c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
			smp.scatter(ra_rcat, dec_rcat, marker=',', c=c, s=5)
			smp.plot(ra, dec, '*', markersize=10, color='orange')
		smp.plot_earth(edges_[0],self.geometry)
		smp.plot_detector(edges_[0],self.geometry,good_detector_list=namelist)
		smp.plot_galactic_plane()
		smp.plot_sum(edges_[0],self.geometry)
		smp.plot_moon(edges_[0],self.geometry)
		smp.plot_continue_source()

		smp.add_subplot(2,1,2)
		smp.title('Energy band 50-300 kev, trigger time ' + utc.fits)
		if self.location50_300[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
			d_chi2 = chi2-chi2.min()
			smp.tricontour(ra_rcat,dec_rcat,d_chi2,[0,2.3,4.61,9.21])
			c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
			smp.scatter(ra_rcat, dec_rcat, marker=',', c=c, s=5)
			smp.plot(ra, dec, '*', markersize=10, color='orange')
		smp.plot_earth(edges_[0],self.geometry)
		smp.plot_detector(edges_[0],self.geometry,good_detector_list=namelist)
		smp.plot_galactic_plane()
		smp.plot_sum(edges_[0],self.geometry)
		smp.plot_moon(edges_[0],self.geometry)
		smp.plot_continue_source()

		smp.savefig(savename)
		smp.close()


class Target(object):

	def __init__(self,detector_names,return_lc_t,candidate,edges_list,new_namelist,location5_50,location50_300,target,geometry,name = None):

		self.detector_names = detector_names
		self.return_lc_t = return_lc_t
		self.candidate = candidate
		self.edges_list = edges_list
		self.new_namelist = new_namelist
		self.location5_50 = location5_50
		self.location50_300 = location50_300
		self.target = target
		self.name = name
		self.geometry = geometry
		self.clock = self.geometry.clock

	def save_triggers(self,savename):

		trig_time = []
		for i in self.edges_list:
			trig_time.append(i[0])
		trig_time = np.array(trig_time)
		utc_time = self.clock.met_to_utc(trig_time).fits
		c = {'met':trig_time,'utc':utc_time,'trig_detector':self.new_namelist}
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename,index=False)

	def save_location_5_50_chi2(self,index,savename):

		if self.location5_50[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
			c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def save_location_5_50(self,index,savename):

		if self.location5_50[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location5_50[index]
			c = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def save_location_50_300_chi2(self,index,savename):

		if self.location50_300[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
			c = {'ra':ra_rcat,'dec':dec_rcat,'chi2':chi2}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def save_location_50_300(self,index,savename):

		if self.location50_300[index] is not None:
			ra,dec,err_r,ra_rcat,dec_rcat,chi2 = self.location50_300[index]
			c = {'ra':[ra],'dec':[dec],'err_r':[err_r]}
			save_perspec = pd.DataFrame(c)
			save_perspec.to_csv(savename,index=False)

	def get_candidate_number(self):

		return len(self.return_lc_t)

	def get_location_5_50(self,index):
		return self.location5_50[index]

	def get_location_50_300(self,index):

		return self.location50_300[index]

	def get_lc_time(self,index):

		return self.return_lc_t[index]

	def get_lc_rate_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_5_50 = lc_arr[1]
		return lc_5_50.T

	def get_lc_bs_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_5_50 = bs_arr[1]
		return bs_5_50.T

	def get_SNR_5_50_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_5_50 = SNR_arr[1]
		return SNR_5_50.T

	def save_candidate_5_50_signal(self,index,savename):

		c = {'met':self.return_lc_t[index],'utc':self.clock.met_to_utc(self.return_lc_t[index]).fits}
		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_ = SNR_arr[1]
		bs_ = bs_arr[1]
		lc_ = lc_arr[1]
		for index,ni in enumerate(self.detector_names):
			c[ni+'_lc'] = lc_[:,index]
			c[ni+'_bs'] = bs_[:,index]
			c[ni+'_SNR'] = SNR_[:,index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename,index=False)

	def get_lc_rate_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_50_300 = lc_arr[2]
		return lc_50_300.T

	def get_lc_bs_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_50_300 = bs_arr[2]
		return bs_50_300.T

	def get_SNR_50_300_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_50_300 = SNR_arr[2]
		return SNR_50_300.T

	def save_candidate_50_300_signal(self, index, savename):

		c = {'met': self.return_lc_t[index], 'utc': self.clock.met_to_utc(self.return_lc_t[index]).fits}
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		SNR_ = SNR_arr[2]
		bs_ = bs_arr[2]
		lc_ = lc_arr[2]
		for index, ni in enumerate(self.detector_names):
			c[ni + '_lc'] = lc_[:, index]
			c[ni + '_bs'] = bs_[:, index]
			c[ni + '_SNR'] = SNR_[:, index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename, index=False)


	def get_lc_rate_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		lc_ = lc_arr[0]
		return lc_.T

	def get_lc_bs_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		bs_ = bs_arr[0]
		return bs_.T

	def get_SNR_for_all_detectors(self,index):

		SNR_arr,lc_arr,bs_arr = self.candidate[index]
		SNR_ = SNR_arr[0]
		return SNR_.T

	def save_candidate_signal(self, index, savename):

		c = {'met': self.return_lc_t[index], 'utc': self.clock.met_to_utc(self.return_lc_t[index]).fits}
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		SNR_ = SNR_arr[0]
		bs_ = bs_arr[0]
		lc_ = lc_arr[0]
		for index, ni in enumerate(self.detector_names):
			c[ni + '_lc'] = lc_[:, index]
			c[ni + '_bs'] = bs_[:, index]
			c[ni + '_SNR'] = SNR_[:, index]
		save_perspec = pd.DataFrame(c)
		save_perspec.to_csv(savename, index=False)

	def plot_candidate_50_300_lc(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		lc_arri = lc_arr[2].T
		bs_arri = bs_arr[2].T
		# SNR_arri = SNR_arr[2].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('Energy band 50-300 kev, trig '+utc.fits)
		for j in range(len(lc_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], lc_arri[j], label='lc ' + detec)
			plt.plot(lc_t - lc_t[0], bs_arri[j], label='bs ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_50_300_SNR(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		# lc_arri = lc_arr[2].T
		# bs_arri = bs_arr[2].T
		SNR_arri = SNR_arr[2].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('Energy band 50-300 kev, trig '+utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], SNR_arri[j], label='SNR ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_5_50_lc(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		lc_arri = lc_arr[1].T
		bs_arri = bs_arr[1].T
		# SNR_arri = SNR_arr[1].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('Energy band 5-50 kev, trig '+utc.fits)
		for j in range(len(lc_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], lc_arri[j], label='lc ' + detec)
			plt.plot(lc_t - lc_t[0], bs_arri[j], label='bs ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_5_50_SNR(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		# lc_arri = lc_arr[1].T
		# bs_arri = bs_arr[1].T
		SNR_arri = SNR_arr[1].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('Energy band 5-50 kev, trig '+utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], SNR_arri[j], label='SNR ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_lc(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		lc_arri = lc_arr[0].T
		bs_arri = bs_arr[0].T
		# SNR_arri = SNR_arr[2].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig '+utc.fits)
		for j in range(len(lc_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], lc_arri[j], label='lc ' + detec)
			plt.plot(lc_t - lc_t[0], bs_arri[j], label='bs ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_candidate_SNR(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		SNR_arr, lc_arr, bs_arr = self.candidate[index]
		# lc_arri = lc_arr[0].T
		# bs_arri = bs_arr[0].T
		SNR_arri = SNR_arr[0].T
		lc_t = self.return_lc_t[index]
		fig = plt.figure(figsize=(10, 10))
		utc = self.clock.met_to_utc(edges_[0])
		fig.suptitle('trig '+utc.fits)
		for j in range(len(SNR_arri)):
			detec = self.detector_names[j]
			plt.subplot(6, 2, j + 1)
			plt.plot(lc_t - lc_t[0], SNR_arri[j], label='SNR ' + detec)
			if detec in namelist:
				plt.axvline(x=edges_[0] - lc_t[0], color='r')
				plt.axvline(x=edges_[-1] - lc_t[0], color='g')
			plt.legend()
		plt.savefig(savename)
		plt.close()

	def plot_sky_map(self, index, savename):

		edges_ = self.edges_list[index]
		namelist = self.new_namelist[index]
		smp = Sky_map(figsize=(10, 10))
		smp.add_subplot(2, 1, 1)
		utc = self.clock.met_to_utc(edges_[0])
		smp.title('Energy band 5-50 kev, trigger time '+utc.fits)
		if self.location5_50[index] is not None:
			ra, dec, err_r, ra_rcat, dec_rcat, chi2 = self.location5_50[index]
			d_chi2 = chi2 - chi2.min()
			smp.tricontour(ra_rcat, dec_rcat, d_chi2, [0, 2.3, 4.61, 9.21])
			c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
			smp.scatter(ra_rcat,dec_rcat,marker = ',',c = c,s = 5)
			smp.plot(ra, dec, '*', markersize=10, color='orange')
		smp.plot_earth(edges_[0], self.geometry)
		smp.plot_detector(edges_[0], self.geometry, good_detector_list=namelist)
		smp.plot_galactic_plane()
		smp.add_source(self.target,self.name)
		smp.plot_sum(edges_[0], self.geometry)
		smp.plot_moon(edges_[0], self.geometry)
		smp.plot_continue_source()

		smp.add_subplot(2, 1, 2)
		smp.title('Energy band 50-300 kev, trigger time '+utc.fits)
		if self.location50_300[index] is not None:
			ra, dec, err_r, ra_rcat, dec_rcat, chi2 = self.location50_300[index]
			d_chi2 = chi2 - chi2.min()
			smp.tricontour(ra_rcat, dec_rcat, d_chi2, [0, 2.3, 4.61, 9.21])
			c = (chi2 - chi2.min()) / (chi2.max() - chi2.min())
			smp.scatter(ra_rcat, dec_rcat, marker=',', c=c, s=5)
			smp.plot(ra, dec, '*', markersize=10, color='orange')
		smp.plot_earth(edges_[0], self.geometry)
		smp.plot_detector(edges_[0], self.geometry, good_detector_list=namelist)
		smp.plot_galactic_plane()
		smp.add_source(self.target,self.name)
		smp.plot_sum(edges_[0], self.geometry)
		smp.plot_moon(edges_[0], self.geometry)
		smp.plot_continue_source()

		smp.savefig(savename)
		smp.close()