import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os


savedir = '/home/laojin/my_lat/cartopy/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)
ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
ax.stock_img()
ny_lon, ny_lat = -75, 43
ax.plot(ny_lon, ny_lat,
         color='blue', linewidth=2, marker='o', transform=ccrs.Geodetic()
         )
plt.savefig(savedir+'A_map.png')
plt.close()