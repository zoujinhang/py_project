
import numpy as np
from GECAM_geometry.satellite import sky_map

import Data_analysis.file as myfile
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import os



link1 = '/home/laojin/my_lat/skyp/tru_locat.txt'
link2 = '/home/laojin/my_lat/skyp/hitl_locat.txt'
savedir = '/home/laojin/my_lat/skyp/'

if os.path.exists(savedir) ==False:
	os.makedirs(savedir)

tru_ra,tru_dec = myfile.readcol(link1)

hitl_ra,hitl_dec = myfile.readcol(link2)

tru = SkyCoord(ra = tru_ra,dec = tru_dec,frame = 'icrs',unit = 'deg')
hitl = SkyCoord(ra = hitl_ra,dec = hitl_dec,frame = 'icrs',unit = 'deg')


sep = []

for i in range(len(tru_ra)):

	sep.append(tru[i].separation(hitl[i]).deg)

sep = np.array(sep)
print(sep)

bins = np.linspace(0,30,29)

bin_n = np.histogram(sep,bins = bins)[0]

bin_n = np.concatenate((bin_n[:1],bin_n))

plt.plot(bins,bin_n)
plt.savefig(savedir+'A_pdf.png')
plt.close()
