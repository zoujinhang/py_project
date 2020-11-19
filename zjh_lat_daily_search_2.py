
from daily_search.satellite import Detectors,Geometry,Locate,Sky_map
from Fermi_tool.daily import Database
from Data_analysis import ch_to_energy,TD_bs,TD_baseline
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt

import os

Fermi_dete = Detectors()

print(Fermi_dete('n0'))

#source = SkyCoord(ra=293.729,dec= 21.3864,frame ='icrs',unit='deg')      #you can change these things
#name = 'SGRJ1935'
#name = 'bn180720598'
#source = SkyCoord(ra=0.53,dec= -2.93,frame ='icrs',unit='deg')
source = SkyCoord(ra=54.5068,dec= -26.9467,frame ='icrs',unit='deg')      #you can change these things
name = 'GRB190114C'

savedir = '/home/laojin/my_lat/location/'
if os.path.exists(savedir)==False:
	os.makedirs(savedir)


#topdir = '/home/laojin/daily_data/'
#topdir = '/media/laojin/TOSHIBA_EXT/daily/'
topdir = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'
#timestart = '2020-04-28T00:19:40.00'
#timestop = '2020-04-28T00:19:47.00'

fermi = Database(topdir)
clock = fermi.clock

#timestart = '2018-07-20T14:19:23'
#timestop = '2018-07-20T14:29:37'
#met_start = clock.utc_to_met(timestart)+1
#met_stop = clock.utc_to_met(timestop)-1


#trigtime = 553789304.654332
trigtime = 569192227


met_start = trigtime - 100
met_stop = trigtime + 150

timestart = clock.met_to_utc(met_start)
timestop = clock.met_to_utc(met_stop)
num = -1
#met_start=609725986.924219+1
#met_stop=609725991.924219-1

binszie = 0.5

bins = np.arange(met_start,met_stop,binszie)
bin_c = 0.5*(bins[1:]+bins[:-1])

indexi = np.where((bin_c>=trigtime)&(bin_c<=trigtime+2))[0]
t_v = bin_c[indexi]
detector_list = fermi.detector
pos = fermi.get_poshist_data(timestart,timestop)
data = fermi.get_detector_data(timestart,timestop)

dete_rate = []
dete_bs = []
t_m = t_v[num]


for ni in fermi.detector:

	ni_ch_n = data[ni]['ch_E']
	ch_n = ni_ch_n['CHANNEL'].values
	e1 = ni_ch_n['E_MIN'].values
	e2 = ni_ch_n['E_MAX'].values

	ni_event = data[ni]['events']
	time = ni_event['TIME'].values
	ch = ni_event['PHA'].values

	time,energy = ch_to_energy(time,ch,ch_n,e1,e2)
	index = np.where((energy>=50)&(energy<=300))[0]
	#index = np.where((energy>=5)&(energy<=50))[0]
	time_e = time[index]

	bin_n = np.histogram(time_e,bins = bins)[0]
	rate = bin_n/binszie
	bin_c,cs,bs =TD_baseline(bin_c,rate)
	#cs,bs = TD_bs(bin_c,rate)
	dete_rate.append(rate[indexi[num]])
	dete_bs.append(bs[indexi[num]])
	plt.plot(bin_c,rate,label='lightcurve')
	plt.plot(bin_c,bs,label = 'bs')
	plt.axvline(x=t_v[num],color = 'r')
	plt.legend()
	plt.savefig(savedir + 'Z_'+ni+'_lightcurve.png')
	plt.close()

print('rate\n',dete_rate)
print('bs \n',dete_bs)



Fermi = Geometry(pos)
print('time_band',Fermi.met_time_band)
print('get_pos',Fermi.get_pos(t_m))
print('get_earth_position',Fermi.get_earth_point(t_m))
print('qsj',Fermi.get_qsj(t_m))
print('detect',Fermi(t_m))
locate = Locate(Fermi)
ra,dec,err_r,ra_rcat,dec_rcat,chi2 = locate.condidate(t_m,np.array(dete_rate),np.array(dete_bs),'case3',[50,300])

print('chi2 max',chi2.max())
print('chi2 min',chi2.min())
print('location',ra,dec,err_r)



c = (chi2-chi2.min())/(chi2.max()-chi2.min())
index_x = np.where(c<1)[0]
smp = Sky_map(figsize=(10,5))
smp.add_subplot(1,1,1)
smp.scatter(ra_rcat[index_x],dec_rcat[index_x],marker = ',',c = c[index_x])
smp.plot_earth(t_v[num],Fermi)
smp.plot_detector(t_v[num],Fermi)
smp.plot(ra,dec,'*',markersize = 10,color = 'orange')
smp.add_source(source,name)
smp.savefig(savedir+'A_sky_map_n_101.png')
smp.close()
