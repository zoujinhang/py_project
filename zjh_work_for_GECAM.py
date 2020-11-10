from zjh_universal_geometry import *
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import os
from mpl_toolkits.basemap import Basemap

def GECAM_geometry_plot(file0,file1,time_edges,savedir,plot_num = 3,points = None,source = None,radius = 10):
	if(os.path.exists(savedir) == False):
		os.makedirs(savedir)
	GECAM_az = np.array([-90.,-50.,-50.,-50.,-50.,-50.,-50.,-16.5,-16.5,-16.5,-16.5,-16.5,-16.5,-16.5,-16.5,-16.5,-16.5,-0.,-0.,-0.,-0.,-0.,-0.,-0.,-0.])*u.deg
	GECAM_zen =np.array([0.,117.83,171.33,224.83,297.83,351.33,44.83,99.5,135.5,171.5,207.5,243.5,279.5,315.5,351.5,27.5,63.5,90.,145.,180.,235.,270.,325.,0.,55.])*u.deg
	name_list = np.arange(1,26,1)
	adq1 = [2,6,12,16,23]
	adq2 = [3,5,11,19,24]
	adq3 = [4,8,15,22,25]
	adq4 = [7,9,13,17,20]
	adq5 = [1,10,14,18,21]
	color1 = []
	for i in name_list:
		if i in adq1:
			color1.append('#2a5caa')#瑠璃色
		elif i in adq2:
			color1.append('#f47920')#橙色
		elif i in adq3:
			color1.append('#1d953f')#薄緑
		elif i in adq4:
			color1.append('#8552a1')#紫
		elif i in adq5:
			color1.append('r')#红
	color2 = []
	for i in name_list:
		if i in adq1:
			color2.append('#f15a22')#金赤
		elif i in adq2:
			color2.append('#aa2116')#绯色
		elif i in adq3:
			color2.append('#009ad6')#青
		elif i in adq4:
			color2.append('#69541b')#国防色
		elif i in adq5:
			color2.append('#843900')#褐色

	fit0 = fits.open(file0)
	data0 = fit0[1].data
	time_0 = data0.field(0)
	fit1 = fits.open(file1)
	data1 = fit1[1].data
	time_1 = data1.field(0)
	time_start = Time_transition().utc_to_met(time_edges[0])
	time_stop = Time_transition().utc_to_met(time_edges[1])
	time_array = np.linspace(time_start, time_stop, plot_num)
	index0 = get_right_index(time_array, time_0)
	index1 = get_right_index(time_array, time_1)

	time_0=time_0[index0]
	qsj1_0= data0.field(1)[index0]
	qsj2_0= data0.field(2)[index0]
	qsj3_0= data0.field(3)[index0]
	qsj4_0= data0.field(4)[index0]
	pos_x_0 = data0.field(8)[index0]
	pos_y_0 = data0.field(9)[index0]
	pos_z_0 = data0.field(10)[index0]

	time_1=time_1[index1]
	qsj1_1= data1.field(1)[index1]
	qsj2_1= data1.field(2)[index1]
	qsj3_1= data1.field(3)[index1]
	qsj4_1= data1.field(4)[index1]
	pos_x_1 = data1.field(8)[index1]
	pos_y_1 = data1.field(9)[index1]
	pos_z_1 = data1.field(10)[index1]

	for i in range(plot_num):
		print('plot ',i)
		time_str = Time_transition().met_to_utc(time_array[i]).iso
		t_0 = time_0[i]
		pos_0 = [pos_x_0[i],pos_y_0[i],pos_z_0[i]]
		qsj_0 = [qsj1_0[i],qsj2_0[i],qsj3_0[i],qsj4_0[i]]
		my_ug0 = U_geometry(qsj_0, GECAM_az, GECAM_zen, sc_pos=pos_0 * u.m, time=t_0,detector_name=name_list)
		t_1 = time_1[i]
		pos_1 = [pos_x_1[i],pos_y_1[i],pos_z_1[i]]
		qsj_1 = [qsj1_1[i],qsj2_1[i],qsj3_1[i],qsj4_1[i]]
		my_ug1 = U_geometry(qsj_1, GECAM_az, GECAM_zen, sc_pos=pos_1 * u.m, time=t_1,detector_name=name_list)
		fig = plt.figure(figsize=(20,25))
		ax0 = fig.add_subplot(2,1,1)
		plt.title('GC1 '+time_str,size = 20)
		mp1 = Basemap(projection = 'moll',lat_0 = 0,lon_0 = 180,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax0)
		my_ug0.detector_plot(radius=radius,points=points,source=source,map = mp1,show_bodies=True,color=color1,detector_style = 'loop')
		ax1 = fig.add_subplot(2,1,2)
		plt.title('GC2 '+time_str,size = 20)
		mp2 = Basemap(projection = 'moll',lat_0 = 0,lon_0 = 180,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax1)
		my_ug1.detector_plot(radius=radius,points=points,source=source,map = mp2,show_bodies=True,color = color1,detector_style = 'cover')
		fig.savefig(savedir+'Z_'+str(i)+'_plot.png')
		plt.close()
		fig = plt.figure(figsize=(20,10))
		ax = fig.add_subplot(1,1,1)
		mp = Basemap(projection = 'moll',lat_0 = 0,lon_0 = 180,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax)
		my_ug0.detector_plot(radius=radius,points=points,source=source,map = mp,show_bodies=False,color=color1,map_style=False,detector_style = 'loop')
		my_ug1.detector_plot(radius=radius,points=points,map = mp,show_bodies=False,color=color1,detector_style = 'cover')
		fig.savefig(savedir+'Z_'+str(i)+'_all_plot.png')
		plt.close()
	return True
def get_right_index(time_array,time_input):
	return_list = []
	for i in time_array:
		t = (time_input-i)**2
		index = np.argmin(t)
		return_list.append(index)
	return np.array(return_list)

topdir = '/home/laojin/shiyan/GECAM_geometry/'
filename0 = 'GECAM01_poshist_all_200702.fits'
filename1 = 'GECAM02_poshist_all_200702.fits'
time = ['2020-07-02T00:00:00','2020-07-02T00:00:05']
link1 = topdir+filename0
link2 = topdir+filename1
GECAM_geometry_plot(link1,link2,time,topdir,plot_num = 3)




