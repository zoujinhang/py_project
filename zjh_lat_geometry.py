
from Data_analysis.geometry import Geometry,Detectors
from Data_analysis import Time_transform
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
datalink = '/home/laojin/trigdata/bn190530430/glg_trigdat_all_bn190530430_v02.fit'
savedir = '/home/laojin/my_lat/geometry/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

hl = fits.open(datalink)
trigtime = hl[0].header['TRIGTIME']
time1 = hl[5].data.field(0)
time2 = hl[5].data.field(1)
q4 = hl[5].data.field(2)
sic = hl[5].data.field(3)
sourec = SkyCoord(116.92,33.95,frame = 'icrs',unit = 'deg')
t_c = (time2+time1)*0.5
dt = 1/24
detectors = Detectors()
fermi_time = Time_transform()
fermi_gbm = Geometry(detector=detectors,time_base=fermi_time)
fermi_gbm.input_pose(quaternion=q4,sc_pos=sic*u.km,time = t_c)
tab = fermi_gbm.get_separation(index = 122,source = sourec)
my_map = fermi_gbm.detector_plot(show_bodies=True,style = 'A',index = 122)
plt.savefig(savedir + 'B_sky_map.png')
plt.close()

n_tab = tab[tab['Detector_index']<=11]
sort_index = np.argsort(n_tab['Separation'])
print(tab.sort('Separation'))
print(n_tab[sort_index]['Detector_index'])
#print(detectors.name_list[(tab[tab['Detector_index']<=11].sort('Separation'))['Detector_index'][:3]])

#fermi_gbm.detector_video(dt,savedir + 'A_sky_map.mp4',show_bodies=True)










