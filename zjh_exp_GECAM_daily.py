
from GECAM_tool.database import Daily_database
from GECAM_tool.satellite import Geometry,Sky_map
from GECAM_tool import Clock
import astropy.units as u
import os
from astropy.coordinates import SkyCoord


#-------------------------------------------------------------------------------
time_list_you_want = ['2020-12-09T21:15:00','2020-12-09T21:16:00',
		      '2020-12-09T21:17:00','2020-12-09T21:18:00']

source = SkyCoord(ra=293.729,dec= 21.3864,frame ='icrs',unit='deg')
name = 'SGRJ1935'

savedir = '/home/laojin/my_lat/GECAM_pipline/'

daily_top = '/home/laojin/GECAM/daily/'	                       #yours

#-------------------------------------------------------------------------------
daily = Daily_database(daily_top)                              #initialize daily database
GC_clock = Clock()                                             #initialize clock transform system

if os.path.exists(savedir)==False:
	os.makedirs(savedir)

met = GC_clock.utc_to_met(time_list_you_want)
utc_t = GC_clock.met_to_utc(met).fits
met_start = met.min()-1
met_stop = met.max()+1

daily_start = GC_clock.met_to_utc(met_start)
daily_stop = GC_clock.met_to_utc(met_stop)

position_data = daily.get_satellite_position(daily_start,daily_stop)

GC_A = Geometry(position_data['A'],pos_unit = u.m, name = 'a',clock=GC_clock)	#initialize GECAM geometry
GC_B = Geometry(position_data['B'],pos_unit = u.m, name = 'b',clock=GC_clock)	#initialize GECAM geometry

A_sep = GC_A.get_separation_with_time(met,source)
A_sep.to_csv(savedir+'A_GECAM_A_'+name+'_sep.csv',index=False)
B_sep = GC_B.get_separation_with_time(met,source)
B_sep.to_csv(savedir+'A_GECAM_B_'+name+'_sep.csv',index=False)

GC_A.get_separation_with_time(met[0],source)
for i,t in enumerate(met):

	smp = Sky_map(figsize = (10,10))
	smp.add_subplot(2,1,1)
	smp.ax.set_title(utc_t[i])
	smp.plot_earth(t,GC_A)
	A_good_dete = GC_A.get_good_detector(t,source)
	smp.plot_continue_source()
	smp.plot_sum(t,GC_A)
	smp.add_source(source,name)
	smp.plot_moon(t,GC_A)
	smp.plot_detector(t,GC_A,good_detector_list=A_good_dete)
	smp.plot_galactic_plane()

	smp.add_subplot(2,1,2)
	smp.plot_continue_source()
	smp.add_source(source,name)
	B_good_dete = GC_B.get_good_detector(t,source)
	smp.plot_earth(t,GC_B)
	smp.plot_sum(t,GC_B)
	smp.plot_moon(t,GC_B)
	smp.plot_detector(t,GC_B,good_detector_list=B_good_dete)
	smp.plot_galactic_plane()
	smp.savefig(savedir + 'B_map_'+str(i)+'.png')
	smp.close()
