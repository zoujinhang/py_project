
from daily_search.satellite import Sky_map
from daily_search import Database
from daily_search import Search
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.coordinates import SkyCoord

datatop = '/media/laojin/TOSHIBA_EXT/daily/'
savetop = '/home/laojin/my_lat/daily_search/new_code1/'
if os.path.exists(savetop) ==False:
	os.makedirs(savetop)

t_start = '2020-04-28T00:00:00'
t_stop = '2020-04-28T01:00:00'
source = SkyCoord(ra=293.729,dec= 21.3864,frame ='icrs',unit='deg')      #you can change these things
name = 'SGRJ1935'
fermi_data = Database(datatop)

data = fermi_data.get_detector_data(t_start,t_stop)
#print('data:',data['n0']['events'])
pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)

search = Search(data,pd_pos_data)
location5_50,location50_300 = search.locate()
#print()
#print(search.candidate)

for index,lc_t in enumerate(search.time_list):

	edges_ = search.edges_list[index]
	namelist = search.name_list[index]
	index_ = np.where((lc_t>=edges_[0]-1)&(lc_t<=edges_[-1]+1))[0]
	t_c = 0.5*(edges_[0]+edges_[-1])
	SNR_arr,lc_arr,bs_arr = search.candidate[index]

	if location5_50[index] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location5_50[index]
		c = (chi2-chi2.min())/(chi2.max()-chi2.min())
		index_x = np.where(c<1)[0]
		smp = Sky_map(figsize=(10,5))
		smp.add_subplot(1,1,1)
		smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = c[index_x],s = 5)
		#smp.tricontour(ra_rcat[index_x],dec_rcat[index_x],chi2[index_x]-chi2.min(),[0,2.3,4.61,9.21])
		smp.tricontour(ra_rcat[index_x], dec_rcat[index_x], chi2[index_x] - chi2.min(), [0,9.21, 18.4, 36.87])
		smp.plot_earth(t_c,search.geometry)
		smp.plot_detector(t_c,search.geometry,good_detector_list=namelist)
		smp.plot(ra,dec,'*',markersize = 10,color = 'orange')
		smp.add_source(source,name)
		smp.savefig(savetop+'A_t'+str(index) + '_sky_map_5_50.png')
		smp.close()

	if location50_300[index] is not None:
		ra,dec,err_r,ra_rcat,dec_rcat,chi2 = location50_300[index]
		c = (chi2-chi2.min())/(chi2.max()-chi2.min())
		index_x = np.where(c<1)[0]
		smp = Sky_map(figsize=(10,5))
		smp.add_subplot(1,1,1)
		smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = c[index_x],s = 5)
		#smp.tricontour(ra_rcat[index_x],dec_rcat[index_x],chi2[index_x]-chi2.min(),[0,2.3,4.61,9.21])
		smp.tricontour(ra_rcat[index_x], dec_rcat[index_x], chi2[index_x] - chi2.min(), [0,9.21, 18.4, 36.87])
		smp.plot_earth(t_c,search.geometry)
		smp.plot_detector(t_c,search.geometry,good_detector_list=namelist)
		smp.plot(ra,dec,'*',markersize = 10,color = 'orange')
		smp.add_source(source,name)
		smp.savefig(savetop+'A_t'+str(index) + '_sky_map_50_300.png')
		smp.close()


	#print(lc_arr)
	for i in range(len(search.energy_band)):
		lc_arri = lc_arr[i].T
		bs_arri = bs_arr[i].T
		SNR_arri = SNR_arr[i].T
		#print(lc_arri)
		plt.figure(figsize = (10,10))

		for j in range(len(lc_arri)):
			detec = search.name[j]
			plt.subplot(6,2,j+1)
			plt.title(detec)
			plt.plot(lc_t[index_]-lc_t[0],lc_arri[j][index_],label = 'lc')
			plt.plot(lc_t[index_]-lc_t[0],bs_arri[j][index_],label = 'bs')
			if detec in namelist:
				plt.axvline(x = edges_[0]-lc_t[0],color = 'r')
				plt.axvline(x = edges_[-1]-lc_t[0],color = 'g')
			plt.legend()

		plt.savefig(savetop + 'A_E'+ str(i) + '_t'+str(index) + '_lc.png')
		plt.close()

		plt.figure(figsize = (10,10))

		for j in range(len(lc_arri)):
			detec = search.name[j]
			plt.subplot(6,2,j+1)
			plt.title(detec)
			plt.plot(lc_t[index_]-lc_t[0],SNR_arri[j][index_],label = 'lc')
			if detec in namelist:
				plt.axvline(x = edges_[0]-lc_t[0],color = 'r')
				plt.axvline(x = edges_[-1]-lc_t[0],color = 'g')
			plt.legend()

		plt.savefig(savetop + 'A_E'+ str(i) + '_t'+str(index) + '_SNR.png')
		plt.close()



