import numpy as np
from daily_search.satellite.detectors import Detectors
from daily_search.perception.triger import check_event,get_subsection_index,get_bins
from daily_search import Database
from daily_search.background import TD_baseline
from daily_search.perception import Event
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
li = ['n3','n6','n7','n8','n9']
#datatop = '/media/laojin/TOSHIBA_EXT/daily/'
datatop = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'
t_start = '2020-12-27T09:00:00.000000000'
t_stop = '2020-12-27T12:00:00.000000000'

fermi_data = Database(datatop)
data = fermi_data.get_detector_data(t_start,t_stop)
clock = fermi_data.clock
wind_list = get_bins(data,wt = 0.1)

met_start = clock.utc_to_met(t_start)

met_stop = clock.utc_to_met(t_stop)

binsize_lightcurve = 0.05
savedir = '/home/laojin/my_lat/location/20020/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

trig = clock.utc_to_met('2020-12-21T05:02:32.8')

for ni in fermi_data.detector:
	plt.figure()
	plt.subplot(1, 1, 1)
	for (wind_start, wind_stop) in wind_list:
		bins_lightcurve = np.arange(wind_start, wind_stop, binsize_lightcurve)
		ni_ch_n = data[ni]['ch_E']
		ni_event = data[ni]['events']
		time = ni_event['TIME'].values
		ch = ni_event['PHA'].values
		lightcurve = Event(ni_ch_n, time, ch)
		t_cc, rate_c = lightcurve(bins=bins_lightcurve, energy_band=[8, 105], channel=True)
		bs = TD_baseline(t_cc, rate_c)
		plt.plot(t_cc, rate_c, label='lightcurve')
		plt.plot(t_cc, bs, label='bs')
	#plt.xlim(trig - 10, trig + 50)
	plt.legend()
	plt.savefig(savedir + 'Z_' + ni + '_lightcurve.png')
	plt.close()
