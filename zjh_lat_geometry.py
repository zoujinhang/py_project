
from Data_analysis.geometry import Geometry,Detectors
from Data_analysis import Time_transform
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u

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

t_c = (time2+time1)*0.5
dt = 1/24
detectors = Detectors()
fermi_time = Time_transform()
fermi_gbm = Geometry(detector=detectors,time_base=fermi_time)
fermi_gbm.input_pose(quaternion=q4,sc_pos=sic*u.km,time = t_c)
my_map = fermi_gbm.detector_plot(show_bodies=True,style = 'A',index = 122)
plt.savefig(savedir + 'B_sky_map.png')
plt.close()
fermi_gbm.detector_video(dt,savedir + 'A_sky_map.mp4',show_bodies=True)










