
import numpy as np
import cartopy.crs as ccrs
import cartopy.geodesic as ccrsge

import matplotlib.pyplot as plt
import os

from copy import copy
import shapely.geometry as sgeom
from spherical_geometry.polygon import SphericalPolygon
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
from cartopy.feature.nightshade import Nightshade
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import datetime
import shapely

def find_side(ls, side):
	"""
	Given a shapely LineString which is assumed to be rectangular, return the
	line corresponding to a given side of the rectangle.
	"""
	minx, miny, maxx, maxy = ls.bounds
	points = {'left': [(minx, miny), (minx, maxy)],
		'right': [(maxx, miny), (maxx, maxy)],
		'bottom': [(minx, miny), (maxx, miny)],
		'top': [(minx, maxy), (maxx, maxy)],}
	return sgeom.LineString(points[side])

def lambert_xticks(ax, ticks):
	"""Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
	te = lambda xy: xy[0]
	lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
	xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
	ax.xaxis.tick_bottom()
	ax.set_xticks(xticks)
	ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])

def lambert_yticks(ax, ticks):
	"""Draw ticks on the left y-axis of a Lamber Conformal projection."""
	te = lambda xy: xy[1]
	lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
	yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
	ax.yaxis.tick_left()
	ax.set_yticks(yticks)
	ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])
	
def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
	"""Get the tick locations and labels for an axis of a Lambert Conformal projection."""
	outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
	axis = find_side(outline_patch, tick_location)
	n_steps = 30
	extent = ax.get_extent(ccrs.PlateCarree())
	_ticks = []
	for t in ticks:
		xy = line_constructor(t, n_steps, extent)
		proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
		xyt = proj_xyz[..., :2]
		ls = sgeom.LineString(xyt.tolist())
		locs = axis.intersection(ls)
		if not locs:
			tick = [None]
		else:
			tick = tick_extractor(locs.xy)
		_ticks.append(tick[0])
	ticklabels = copy(ticks)
	while True:
		try:
			index = _ticks.index(None)
		except ValueError:
			break
		_ticks.pop(index)
		ticklabels.pop(index)
	return _ticks, ticklabels
	
	
	

savedir = '/home/laojin/my_lat/cartopy/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)






plt.figure(figsize=(20,10))
#my_globe = ccrs.Globe(inverse_flattening=-1)
#ax = plt.subplot(111,projection=ccrs.Mollweide(central_longitude=180))
ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180),facecolor = '#f6f5ec')#第二个参数修改背景色

#ax.stock_img()
#ax.coastlines(resolution='110m')

lats_y = np.arange(-75,76,15)
lats_x = np.zeros(lats_y.size)

print(ax.get_extent())

#ax.set_yticks(lats_y,crs=ccrs.Mollweide(central_longitude=180))
#ax.gridlines()
#ax.spines['bottom'].set_position(('data',0))
#------------------------------------------------
#ax.set_xlim(ax.get_extent()[:2])
#ax.set_ylim(ax.get_extent()[2:])
ax.set_global()                                            #跟上面两句效果一样
#------------------------------------------------
ax.invert_xaxis()                                          #设置轴反向！！牛逼！！，这个轴翻转要放到最后
#ax.set_ylim([-90,90])


ny_lon, ny_lat = -75, 43
ax.plot(ny_lon, ny_lat,color='blue', linewidth=2, marker='o', transform=ccrs.Geodetic())
ax.plot(0,0,color = 'r',transform=ccrs.Geodetic(), linewidth=2, marker='o')
ax.plot(180,0,color = 'g', linewidth=2, marker='o',transform=ccrs.Geodetic())
ax.plot(270,0,color = 'k', linewidth=2, marker='o',transform=ccrs.Geodetic())
ax.plot([260,270,280,290,330,350,0,10],[-60,-50,-40,-30,-20,0,10,20],'-',color = 'k', linewidth=2,transform=ccrs.Geodetic())
#ax.plot(x,y,'-',color = 'r',linewidth = 2,transform = ccrs.Geodetic())

poly = SphericalPolygon.from_cone(0,0,60,steps=100)
print(np.array(list(poly.to_radec())[0]).T)
x,y = [p for p in poly.to_radec()][0]
new_xy = list(zip(x,y))
print('new_xy:\n',new_xy)
new_xy2 = np.array(list(poly.to_radec())[0]).T
circle_points = ccrsge.Geodesic().circle(lon=0, lat=0, radius=60, n_samples=200, endpoint=False)
print('list:',list(circle_points[0]))
geom = shapely.geometry.Polygon(new_xy2)
print('array:',np.array(list(circle_points)))
print(geom)
ax.add_geometries((geom,), crs=ccrs.Mollweide(central_longitude=180), facecolor='red', edgecolor='none', linewidth=0)
#poly1 = Polygon(new_xy[::-1],facecolor='#74787c',edgecolor='#74787c',linewidth=2, alpha=0.5,transform=ccrs.Geodetic())
poly1 = Polygon(new_xy2[::-1],facecolor='#74787c',edgecolor='#74787c',linewidth=2, alpha=0.5,transform=ccrs.Geodetic())
poly2 = Polygon(np.array(list(circle_points)),facecolor='#74787c',edgecolor='#74787c',linewidth=2, alpha=0.5,transform=ccrs.Geodetic())
ax.add_patch(poly2)
#ax.add_geometries(poly1)
#ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
#ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
#ax.add_patch(poly1)
#ax.add_patch(mpatches.Circle(xy=[0, 0], radius=60, color='red', alpha=0.3, transform=ccrs.Geodetic(), zorder=0,joinstyle='bevel'))
#ax.add_patch(mpatches.CirclePolygon(xy=[0, 0], radius=60, resolution=100, color='g', alpha=0.3, transform=ccrs.Geodetic(), zorder=0,joinstyle='miter'))
#lambert_xticks(ax, xticks)
#lambert_yticks(ax, yticks)

lons_x = np.arange(0,360,30)
lons_y = np.zeros(lons_x.size)

date = datetime.datetime(1999, 12, 31, 12)
ax.add_feature(Nightshade(date, alpha=0.2))
#print(Nightshade(date, alpha=0.2))
xticks = list(range(-180, 180, 30))
yticks = list(range(-90, 90, 15))
ax.gridlines(xlocs=xticks, ylocs=yticks)

y = np.linspace(-90,90,100)
y_z = np.zeros(y.size)

x = np.linspace(0,360,100)
x_z = np.zeros(x.size)
ax.text(10000000,10000000,'yayayayaya',size = 20)

proj_xyz = ax.projection.transform_points(ccrs.Geodetic(),np.array([0,30]), np.array([0,0]))
ax.text(proj_xyz[0][0],proj_xyz[0][1],'240+15',size = 20,ha = 'right')
dx_ = np.abs(proj_xyz[0][0]-proj_xyz[1][0])
#print(proj_xyz)

lats_y_ticke = ax.projection.transform_points(ccrs.Geodetic(),lats_x+180+179.99, lats_y*1.1)
lats_y_x = lats_y_ticke[:,0]*0.86
lats_y_y = lats_y_ticke[:,1]
#print(lats_y_ticke)
#print(lats_y_x)
#print(lats_y_y)
for index,i in enumerate(lons_x):
	
	#ax.plot(y_z+i,y,'-',color = '#d3d7d4',linewidth = 1,transform = ccrs.Geodetic())
	ax.text(i,lons_y[index],r'$%d^{\circ}$'%i,transform = ccrs.Geodetic(),size = 20)

for index,i in enumerate(lats_y):
	
	#ax.plot(x,x_z+i,'-',color = '#d3d7d4',linewidth = 1,transform = ccrs.Geodetic())
	ax.text(lats_y_x[index]+dx_,lats_y_y[index],r'$%d^{\circ}$'%i,size = 20,ha = 'right',va = 'center')
	
ax.set_xlabel('sssssss')
plt.savefig(savedir+'A_map.png')
plt.close()




plt.figure()
plt.subplot(111,projection = "mollweide")
plt.title('Mollweide')
plt.grid(True)
plt.savefig(savedir + 'B_map.png')
plt.close()





