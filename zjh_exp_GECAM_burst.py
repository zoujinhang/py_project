
from GECAM_tool.database import Burst_database
from GECAM_tool.satellite import Geometry,Sky_map
from GECAM_tool import Clock
import astropy.units as u
import os
from astropy.coordinates import SkyCoord
import numpy as np

#-------------------------------------------------------------------------

source = SkyCoord(ra=293.729,dec= 21.3864,frame ='icrs',unit='deg')      #you can change these things
name = 'SGRJ1935'
samplename = 'bn201209_213202'

savedir = '/home/laojin/my_lat/GECAM_pipline/'

topdir = '/home/laojin/GECAM/bursts/'

#-------------------------------------------------------------------------

if os.path.exists(savedir)==False:
	os.makedirs(savedir)

GC_clock = Clock()



catalog = Burst_database(topdir)
position_data = catalog.get_satellite_position(samplename)
GC_A = Geometry(position_data['A'],pos_unit = u.m, name = 'a',clock=GC_clock)


#--------------------------------------------------------------------------
met = np.linspace(GC_A.met_time_band[0],GC_A.met_time_band[1],10)       #you can change these things
#--------------------------------------------------------------------------


utc_t = GC_clock.met_to_utc(met).fits
A_seq = GC_A.get_separation_with_time(met,source)
A_seq.to_csv(savedir+'A_GECAM_A_'+name+'_sep.csv',index=False)
GC_B = Geometry(position_data['B'],pos_unit = u.m, name = 'b',clock=GC_clock)
B_seq = GC_B.get_separation_with_time(met,source)
B_seq.to_csv(savedir+'A_GECAM_B_'+name+'_sep.csv',index=False)
for i in range(10):

	t = GC_A.met_time_band[0]+10*i+1
	utc = GC_clock.met_to_utc(t).fits
	smp = Sky_map(figsize = (10,10))
	smp.add_subplot(2,1,1)
	smp.ax.set_title(utc)
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
