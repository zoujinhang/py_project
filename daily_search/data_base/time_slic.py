
from astropy.time import Time
import pandas as pd

def time_slic(time_start,time_stop,format = None,scale='utc',H = 3):

	if isinstance(time_start,Time) ==False:
		time_start = Time(time_start, format=format, scale=scale)
	if isinstance(time_stop,Time) ==False:
		time_stop =Time(time_stop,format=format,scale=scale)
	re = pd.date_range(time_start.fits,time_stop.fits,freq = str(H)+'H')
	re = pd.DataFrame({'time':re})
	stp = pd.date_range(time_stop.fits,freq = 'H',periods=1)
	stp = pd.DataFrame({'time':stp})

	t_ar = pd.concat([re,stp])
	t_ar.drop_duplicates('time','first',inplace=True,ignore_index=True)
	t_ar.sort_values(by = 'time',inplace=True,ignore_index=True)
	time_list = t_ar['time'].values
	t_ar1 = pd.DataFrame({'time_start':time_list[:-1],'time_stop':time_list[1:]})
	return t_ar1.values
