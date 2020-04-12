import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from multiprocessing import Pool
import Data_analysis.file as myfile
from Data_analysis.geometry import Geometry,Detectors
from Data_analysis import Time_transform,Separate_source,ch_to_energy,TD_baseline
from Data_analysis import get_bayesian_flash,get_bayesian_duration,background_correction
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import bayesian_blocks


databaselink = '/media/laojin/Elements/trigdata/'
yearlist = [2016]


def get_sample_dir_list(yearlist,databaselink):
	sample_dir_list = []
	for year in yearlist:
		topdir = databaselink + str(year) +'/'
		dirlist1 = os.listdir(topdir)
		for dirl in dirlist1:
			if os.path.isdir(topdir + dirl):
				sample_dir_list.append([topdir + dirl+'/',dirl])
	return sample_dir_list



def analysis_one_sample(input_list):
	'''
	
	:param input_list:
	:return:
	'''
	
	NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
	BGO = ['b0','b1']
	files = get_file(input_list,NaI,BGO)
	
	if (files['trigdat'] is not None) and (files['loc'] is not None):
		trigtime = files['trigdat'][0].header['TRIGTIME']
		fermi_gbm = get_fermi_geometry(files['trigdat'])
		detector = fermi_gbm.detectors
		t_c = fermi_gbm.Time_transition(fermi_gbm.time,astropyTime = True)
		loc_hl = files['loc']
		ra = loc_hl[0].header['RA_OBJ']
		dec = loc_hl[0].header['DEC_OBJ']
		source = SkyCoord(ra,dec,frame = 'icrs',unit = 'deg')
		trigtime_index = np.argmin((t_c-trigtime)**2)
		
		tab = fermi_gbm.get_separation(trigtime_index,source=source)      #获得夹角
		ni_tab = tab[tab['Detector_index']<=11]
		sort_index = np.argsort(ni_tab['Separation'])
		good_ni =detector.name_list[ni_tab[sort_index]['Detector_index']] #获得NaI夹角最小探头名列表
		#good_ni = ni_tab.sort('Separation')[:3]#这个排序有问题，问题不明
		bgoi_tab = tab[tab['Detector_index']>11]
		sort_index = np.argsort(bgoi_tab['Separation'])
		good_bgo = detector.name_list[bgoi_tab[sort_index]['Detector_index']] #获得BGO夹角最小探头名列表
		
		
		
	else:
		pass
def light_curve_analysis(file,NaI,BGO,good_ni,good_bi,txtsave,plotsave):
	dt = 0.064
	for ni in NaI:
		if file[ni] is not None:
			hl = file[ni]
			trigtime = hl[0].header['TRIGTIME']
			time = hl[2].data.field(0)
			ch = hl[2].data.field(1)
			t = time - trigtime
			ch_index = np.where((ch>=8)&(ch<110))[0]
			t = t[ch_index]
			#ch= ch[ch_index]
			bins = np.arange(t[0],t[-1],dt)
			bin_n,bin_edges = np.histogram(t,bins = bins)
			t_c = (bin_edges[1:] + bin_edges[:-1])*0.5
			rate = bin_n/dt
			t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
			if ni in good_ni:
				rate_sm = cs_rate+bs_rate.mean()
				bin_n_sm = np.round(rate_sm*dt)
				edges = bayesian_blocks(t_c,bin_n_sm,fitness='events',p0 = 0.001)
				result = background_correction(t_c,rate_sm,edges,degree = 50)
				startedges,stopedges = get_bayesian_duration(result,sigma = 5)
				if len(startedges)>0:
					dir_,file = os.path.split(txtsave[0])
					if os.path.exists(dir_) == False:
						os.makedirs(dir_)
					myfile.printdatatofile(txtsave[0],data = [startedges,stopedges],format = ['.5f','.5f'])
					flash_start,flash_stop = get_bayesian_flash(result,startedges,stopedges)
					dir_,file = os.path.split(txtsave[1])
					if os.path.exists(dir_) == False:
						os.makedirs(dir_)
					myfile.printdatatofile(txtsave[1],data = [flash_start,flash_stop],format = ['.5f','.5f'])
				

def get_fermi_geometry(hl):
	
	time1 = hl[5].data.field(0)
	time2 = hl[5].data.field(1)
	q4 = hl[5].data.field(2)
	sic=hl[5].data.field(3)
	t_c = (time2 + time1)*0.5
	detectors = Detectors()
	fermi_time = Time_transform()
	fermi_gbm=Geometry(detector=detectors,time_base = fermi_time)
	fermi_gbm.input_pose(quaternion=q4,sc_pos=sic*u.km,time = t_c)
	return fermi_gbm

def plot_count_map(file,NaI,BGO,savedirlist):
	
	plt.figure(figsize=(40,40))
	n = 1
	for ni in NaI:
		if file[ni] is not None:
			hl = file[ni]
			trigtime = hl[0].header['TRIGTIME']
			time = hl[2].data.field(0)
			ch = hl[2].data.field(1)
			t = time - trigtime
			ch_n = hl[1].data.field(0)
			e1 = hl[1].data.field(1)
			e2 = hl[1].data.field(2)
			plt.subplot(4, 4, n)
			plt.title(ni)
			ni_event = Separate_source(t,ch,ch_n)
			s_t,s_ch = ni_event.get_S_t_and_ch()
			new_t,new_energy = ch_to_energy(s_t,s_ch,ch_n,e1,e2)
			plt.plot(new_t,new_energy,',',color = 'k')
			plt.xlabel('time (s)')
			plt.ylabel('Energy (KeV)')
			plt.xlim(new_t[0],new_t[-1])
			plt.ylim(8,900)
			plt.yscale('log')
			n = n+1
	for bi in BGO:
		if file[bi] is not None:
			hl = file[bi]
			trigtime = hl[0].header['TRIGTIME']
			time = hl[2].data.field(0)
			ch = hl[2].data.field(1)
			t = time - trigtime
			ch_n = hl[1].data.field(0)
			e1 = hl[1].data.field(1)
			e2 = hl[1].data.field(2)
			plt.subplot(4, 4, n)
			plt.title(bi)
			ni_event = Separate_source(t,ch,ch_n)
			s_t,s_ch = ni_event.get_S_t_and_ch()
			new_t,new_energy = ch_to_energy(s_t,s_ch,ch_n,e1,e2)
			plt.plot(new_t,new_energy,',',color = 'k')
			plt.xlabel('time (s)')
			plt.ylabel('Energy (KeV)')
			plt.xlim(new_t[0],new_t[-1])
			plt.ylim(800,20000)
			plt.yscale('log')
			n = n+1
	for savedir in savedirlist:
		dir_,file = os.path.split(savedir)
		if os.path.exists(dir_) == False:
			os.makedirs(dir_)
		plt.savefig(savedir)
	plt.close()




def plot_sky_map(fermi_gbm,source,trigtime, savedirlist):
	
	t_c = fermi_gbm.Time_transition(fermi_gbm.time,astropyTime = True)
	index_0 = np.argmin((t_c-trigtime)**2)
	fermi_gbm.detector_plot(show_bodies = True,index = index_0,source = source)
	for savedir in savedirlist:
		dir_,file = os.path.split(savedir)
		if os.path.exists(dir_) == False :
			os.makedirs(dir_)
		plt.savefig(savedir)
	plt.close()
	
def get_file(input_list,NaI,BGO):
	sample = input_list[1]
	sampledir = input_list[0]
	data = {}
	loc_name_list = myfile.findfile(sampledir,'glg_locprob_all_'+sample+'_v*')
	if len(loc_name_list) >= 1:
		data['loc'] = fits.open(sampledir+loc_name_list[0])
	else:
		data['loc'] = None
	trigdat_name_list = myfile.findfile(sampledir,'glg_trigdat_all_'+sample+'_v*')
	if len(trigdat_name_list)>=1:
		data['trigdat'] = fits.open(sampledir+trigdat_name_list[0])
	else:
		data['trigdat'] = None
	for ni in NaI:
		ni_list = myfile.findfile(sampledir,'glg_tte_'+ni+'_'+sample+'_v*')
		if len(ni_list)>=1:
			data[ni] = fits.open(sampledir + ni_list[0])
		else:
			data[ni] = None
	for bi in BGO:
		bi_list = myfile.findfile(sampledir,'glg_tte_'+bi+'_'+sample+'_v*')
		if len(bi_list)>=1:
			data[bi] = fits.open(sampledir+bi_list[0])
		else:
			data[bi] = None
	return data
#print(get_sample_dir_list(yearlist,databaselink))


