from Fermi_tool import *
import os
import numpy as np
import matplotlib.pyplot as plt
import Data_analysis.file as myfile




savedir = '/home/laojin/my_lat/spectrum/'
data_top = '/media/laojin/Elements/trigdata/'

if os.path.exists(savedir)  ==False:
	os.makedirs(savedir)

def get_dir_list(year,sample,databaselink):
	sample_dir_list = []
	for i in range(len(year)):
		topdir = databaselink + str(year[i]) +'/' +sample[i]+'/'
		sample_dir_list.append([topdir,sample[i]])
	return sample_dir_list
	
def analysis_one_sample(input_list):
	'''
	
	:param input_list:
	:return:
	'''
	
	NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
	BGO = ['b0','b1']
	sample = input_list[1]
	files = get_file(input_list,NaI,BGO)
	savetop = '/home/laojin/my_lat/spectrum/'
	all_sky_map = savetop + '/A_skymap/D_'+sample+'_skymap.png'
	sky_map = savetop  + '/' + sample + '/D_skymap.png'
	
	all_light_curve_savedir = savetop  + '/A_lightcurve/A_'+sample+'_lightcurve.png'
	light_curve_savedir = savetop  + '/' + sample + '/A_all_lightcurve.png'
	
	txt_savedir = savetop  + '/' + sample +'/'
	
	all_duration_savedir = savetop  + '/A_duration/B_'+sample+'_duration.png'
	duration_savedir = savetop  + '/' + sample + '/B_duration.png'
	
	allcountmapdir = savetop  + '/A_count_map/C_'+sample+'_countmap.png'
	countmapdir = savetop  + '/' + sample + '/C_countmap.png'
	
	if (files['trigdat'] is not None) and (files['loc'] is not None):
		trigtime = files['trigdat'][0].header['TRIGTIME']
		fermi_gbm = get_fermi_geometry(files['trigdat'])#构建gbm的空间几何
		detector = fermi_gbm.detectors#读取探测器
		t_c = fermi_gbm.Time_transition.batch_utc_to_met(fermi_gbm.time,astropyTime = True)
		
		loc_hl = files['loc']
		ra = loc_hl[0].header['RA_OBJ']
		dec = loc_hl[0].header['DEC_OBJ']
		source = SkyCoord(ra,dec,frame = 'icrs',unit = 'deg')
		plot_sky_map(fermi_gbm,source,trigtime,[all_sky_map,sky_map])#--------skymap
		trigtime_index = np.argmin((t_c-trigtime)**2)
		
		tab = fermi_gbm.get_separation(trigtime_index,source=source)      #获得夹角
		ni_tab = tab[tab['Detector_index']<=11]
		sort_index = np.argsort(ni_tab['Separation'])
		good_ni =detector.name_list[ni_tab[sort_index]['Detector_index']] #获得NaI夹角最小探头名列表
		#good_ni = ni_tab.sort('Separation')[:3]#这个排序有问题，问题不明
		bgoi_tab = tab[tab['Detector_index']>11]
		sort_index = np.argsort(bgoi_tab['Separation'])
		good_bgo = detector.name_list[bgoi_tab[sort_index]['Detector_index']] #获得BGO夹角最小探头名列表
		print('plot light curve...')
		light_curve_analysis(files,NaI,BGO,good_ni[:3],good_bgo[:1],txt_savedir,
		                     [all_duration_savedir,duration_savedir],
		                     [all_light_curve_savedir,light_curve_savedir],
		                     txx=True,WT = True)
		print('plot count map...')
		plot_count_map(files,NaI,BGO,[allcountmapdir,countmapdir])
		
	else:
		pass

years,samples = myfile.readcol('/home/laojin/my_lat/spectrum/sample_list.txt')

input_list = get_dir_list(years,samples,data_top)

for value in input_list:
	analysis_one_sample(value)














