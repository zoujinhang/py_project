import pandas as pd
import os
from GECAM_tool.satellite import Sky_map
from astropy.coordinates import SkyCoord

savedir = '/home/laojin/my_work/GRB201021_search/'
datalink = '/home/laojin/my_work/GRB201021_search/browse_results.xls'

if os.path.exists(savedir) == False:

	os.makedirs(savedir)

target = SkyCoord(ra = 128.426,dec = 27.712,frame = 'icrs',unit= 'deg')

df = pd.read_excel(datalink,sheet_name = 'fermigbrst')

grbs_location = SkyCoord(ra = df['ra'].values,dec = df['dec'].values,frame = 'icrs',unit = 'deg')

sep = grbs_location.separation(target).deg

good_index = sep<9

good_grb = df[good_index]
good_grb.to_csv(savedir +'B_good_grbs.csv',index=False)
print(grbs_location)
smp = Sky_map(figsize = (10,5))
smp.plot(df['ra'].values,df['dec'].values,'.')
smp.add_source(target)
smp.plot(good_grb['ra'].values,good_grb['dec'].values,'.',color = 'r')
smp.plot_galactic_plane()
smp.savefig(savedir + 'A_sky_map.png')
smp.close()





