import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import cartesian_to_spherical,SkyCoord,get_sun,get_body_barycentric
from ..Time_transform import Time_transform
from .detectors import Detectors
from spherical_geometry.polygon import SphericalPolygon
from astropy.table import Table
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
from matplotlib.patches import Polygon

class Geometry(object):
	
	def __init__(self,detector=None,time_base = None):
		'''
		
		:param detector:
		:param time_base:
		'''
		if time_base is None:
			self.Time_transition = Time_transform()
		else:
			self.Time_transition = time_base
		if detector is None:
			self.detectors = Detectors()
		else:
			self.detectors = detector
	
		
	def input_pose(self,quaternion,sc_pos = None,time = None,ef_radius = 60):
		'''
		
		:param quaternion:[[q1,q2,q3,q0],[q1,q2,q3,q0]]
		:param sc_pos:[[x,y,z],[x,y,z]]
		:param time:
		:param ef_radius:
		:return:
		'''
		if time is not None:
			if (isinstance(time, str)):
				self.time = Time(time)
			else:
				self.time = self.Time_transition.batch_met_to_utc(time)
		else:
			self.time = None
		self.radius = ef_radius
		self.quaternion = quaternion
		self.index = np.arange(len(self.quaternion),dtype = int)
		self.sc_pos = sc_pos
		self.detectors.input_quaternion(self.quaternion)
		#self.all_sky = self.all_sky(num_points = 2000)
		
	def all_sky(self,num_points = 3000):
		
		ra_grid_tmp = np.linspace(0, 360, num_points)
		dec_range = [-90, 90]
		cosdec_min = np.cos(np.deg2rad(90.0 + dec_range[0]))
		cosdec_max = np.cos(np.deg2rad(90.0 + dec_range[1]))
		v = np.linspace(cosdec_min, cosdec_max, num_points)
		v = np.arccos(v)
		v = np.rad2deg(v)
		v -= 90.
		dec_grid_tmp = v
		ra_grid = np.zeros(num_points ** 2)
		dec_grid = np.zeros(num_points ** 2)
		itr = 0
		for ra in ra_grid_tmp:
			for dec in dec_grid_tmp:
				ra_grid[itr] = ra
				dec_grid[itr] = dec
				itr += 1
		return SkyCoord(ra=ra_grid, dec=dec_grid, frame='icrs', unit='deg')
	
		
	def get_detector_centers(self,index = None):
		if index is not None:
			retur_list = []
			for i in index:
				retur_list.append(self.detectors.center_all[i])
			return retur_list
		else:
			return self.detectors.center_all
	
	def get_detector_index(self,index = None):
		if index is not None:
			retur_list = []
			for i in index:
				retur_list.append(self.detectors.detector_index[i])
			return retur_list
		else:
			return self.detectors.detector_index
	
	def get_good_detector_centers(self,source = None,index = None):
		'''
		
		:param source:
		:return:
		'''
		
		if source is not None:
			good_detector_centers = []
			good_detector_index = []
			if index is not None:
				center_all_list = []
				for i in index:
					center_all_list.append(self.detectors.center_all[i])
				
			else:
				center_all_list = self.detectors.center_all
			for center_all in center_all_list:
				condition = center_all.separation(source)<=self.radius * u.degree
				good_detector_index.append(condition)
				good_detector_centers.append(center_all[condition])
			return good_detector_index, good_detector_centers
		else:
			print('No source! return []')
			return None
	
	
	def get_separation(self,index,source = None):
		'''
		
		:param source:
		:return:
		'''
		tab = Table(names=["Detector_index", "Separation"], dtype=[int, np.float64])
		if source is not None:
			for index1,center in enumerate(self.detectors.center_all[index]):
				sep = center.separation(source)
				tab.add_row([index1,sep])
			tab['Separation'].unit = u.degree
			#tab.sort('Separation')
			return tab
		else:
			print('No source! return None')
			return None

	def get_fov(self,conter,radius = 10.):
		fov_point_list = []
		for conter_i in conter:
			poly = SphericalPolygon.from_cone(conter_i.ra.value,conter_i.dec.value,radius,steps=100)
			x,y = [p for p in poly.to_radec()][0]
			fov_point_list.append([radius,x,y])
		return fov_point_list
		
		
	def get_earth_point(self,index = None):
		if self.sc_pos is not None:
			earth_point_list = []
			if index is not None:
				sc_pos_list = self.sc_pos[index]
			else:
				sc_pos_list = self.sc_pos
			for sc_pos in sc_pos_list:
				position = cartesian_to_spherical(-sc_pos[0],-sc_pos[1],-sc_pos[2])
				xyz_position = SkyCoord(position[2].deg,position[1].deg,frame='icrs',unit='deg')
				earth_radius = 6371. * u.km
				fermi_radius = np.sqrt((sc_pos ** 2).sum())
				radius_deg = np.rad2deg(np.arcsin((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value)
				poly = SphericalPolygon.from_cone(position[2].deg,position[1].deg,radius_deg,steps=100)
				x,y = [p for p in poly.to_radec()][0]
				earth_point_list.append([xyz_position,radius_deg,x,y])
			return earth_point_list
		else:
			print('No satellite position!')
			return None
			
	def detector_plot(self,radius = 10.0,source=None,points = None,good = False,projection = 'moll',
	                  lat_0 = 0,lon_0 = 180,ax = None,show_bodies = False,
	                  style = 'A',index = None):
		pole = SkyCoord([0, 0], [90, -90], frame='icrs', unit='deg')
		if ax is None:
			fig = plt.figure(figsize = (20,10))
			ax = fig.add_subplot(1,1,1)
			ax.set_title(str(index))
		map = Basemap(projection=projection,lat_0=lat_0,lon_0 = lon_0,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax)

		if good and source :
			index_,centor = self.get_good_detector_centers(source,index = [index])
			index_ = index_[0]
			centor = centor[0]
		else:
			index_ = self.get_detector_index(index=[index])[0]
			centor = self.get_detector_centers(index = [index])[0]
		
		if show_bodies and self.sc_pos is not None:
			#print('plot_earth!')
			if projection in ['moll']:
				postion, r, lon, lat = self.get_earth_point(index=[index])[0]
				lon_lis, lat_lis = get_poly(postion, r, lon, lat, pole,lon_0)
				for i in range(len(lon_lis)):
					x, y = map(lon_lis[i], lat_lis[i])
					earth = Polygon(list(zip(x, y)), facecolor='#90d7ec', edgecolor='#90d7ec',
					                linewidth=0, alpha=1)
					ax.add_patch(earth)
			else:
				postion, r, lon, lat = self.get_earth_point(index=[index])[0]
				lon, lat = map(lon, lat)
				earth = Polygon(list(zip(lon, lat)), facecolor='#90d7ec', edgecolor='#90d7ec',
				                linewidth=2, alpha=1)
				ax.add_patch(earth)
			if self.time is not None:
				#print(self.time)
				earth_r = get_body_barycentric('earth', self.time[index])
				moon_r = get_body_barycentric('moon', self.time[index])
				r_e_m = moon_r - earth_r
				r = self.sc_pos[index] - np.array([r_e_m.x.value, r_e_m.y.value, r_e_m.z.value]) * u.km
				moon_point_d = cartesian_to_spherical(-r[0], -r[1], -r[2])
				moon_ra, moon_dec = moon_point_d[2].deg, moon_point_d[1].deg
				moon_point = SkyCoord(moon_ra, moon_dec, frame='icrs', unit='deg')
				moon_ra, moon_dec = map(moon_point.ra.deg, moon_point.dec.deg)
				map.plot(moon_ra, moon_dec, 'o', color='#72777b', markersize=20)
				plt.text(moon_ra, moon_dec - 800000, 'moon', size=20)
			if show_bodies and self.time is not None:
				tmp_sun = get_sun(self.time[index])
				sun_position = SkyCoord(tmp_sun.ra.deg, tmp_sun.dec.deg, unit='deg', frame='icrs')
				sun_ra, sun_dec = map(sun_position.ra.value, sun_position.dec.value)
				map.plot(sun_ra, sun_dec, 'o', color='#ffd400', markersize=40)
				plt.text(sun_ra - 550000, sun_dec - 200000, 'sun', size=20)
		
		fovs = self.get_fov(centor,radius)
		if projection in ['moll']:
			for i,v in enumerate(index_):
				r,ra,dec = fovs[i]
				lon_lis, lat_lis = get_poly(centor[i],r,ra,dec,pole,lon_0)
				#print(str(self.detectors.name_list[v]))
				#print(lon_lis, lat_lis)
				for ij in range(len(lon_lis)):
					x,y = map(lon_lis[ij],lat_lis[ij])
					detec = Polygon(list(zip(x,y)),facecolor=self.detectors.color_list[v],edgecolor=self.detectors.color_list[v],linewidth=2, alpha=0.5)
					ax.add_patch(detec)
				ra_x,dec_y = map(centor[i].ra.value+2.5,centor[i].dec.value-1)
				plt.text(ra_x, dec_y,str(self.detectors.name_list[v]), color=self.detectors.color_list[v], size=22)
		else:
			for i,v in enumerate(index_):
				r,ra,dec = fovs[i]
				detec = Polygon(list(zip(ra,dec)),facecolor=self.detectors.color_list[v],edgecolor=self.detectors.color_list[v],linewidth=2, alpha=0.5)
				ax.add_patch(detec)
				ra_x,dec_y = map(centor[i].ra.value+2.5,centor[i].dec.value-1)
				plt.text(ra_x, dec_y,str(self.detectors.name_list[v]), color=self.detectors.color_list[v], size=22)
		
		if source:
			ra, dec = map(source.ra.value, source.dec.value)
			map.plot(ra, dec, '*', color='#f36c21', markersize=20.)
		if points:
			ra, dec = map(points.ra.value, points.dec.value)
			map.plot(ra, dec, '*', color='#c7a252', markersize=20.)
		
		if projection == 'moll':
			az1 = np.arange(0, 360, 30)
			zen1 = np.zeros(az1.size) + 2
			azname = []
			for i in az1:
				azname.append(r'${\/%s\/^{\circ}}$' % str(i))
			x1, y1 = map(az1, zen1)
			for index1, value in enumerate(az1):
				plt.text(x1[index1], y1[index1], azname[index1], size=20)
		map.drawmeridians(np.arange(0, 360, 30),dashes=[1,0],color='#d9d6c3')
		map.drawparallels(np.arange(-90, 90, 15),dashes=[1,0], labels=[1,0,0,1], color='#d9d6c3',size = 20)
		map.drawmapboundary(fill_color='#f6f5ec')
		return map
		
	def detector_video(self,dt,savevdir,radius = 10.0,source=None,points = None,good = False,projection = 'moll',
	                   lat_0 = 0,lon_0 = 180.,ax = None,show_bodies = False,style = 'A'):
		pole = SkyCoord([0, 0], [90, -90], frame='icrs', unit='deg')
		name_list = self.detectors.name_list
		color_list = self.detectors.color_list
		
		if good and source :
			index_list,centor_list = self.get_good_detector_centers(source)
		else:
			index_list = self.get_detector_index()
			centor_list = self.get_detector_centers()
		if show_bodies and self.sc_pos is not None:
			sc_pos = self.sc_pos
			earth_points_list = self.get_earth_point()
		else:
			sc_pos = None
			earth_points_list = None
		if self.time is not None:
			time = self.time
		else:
			time = None
		#fig,ax = plt.subplots()
		fig = plt.figure(figsize=(20, 10))
		ax = fig.add_subplot(111)
		duration = dt * len(self.index)
		
		def make_frame(t):
			n = int(t/dt)
			ax.clear()
			ax.set_title(str(n))
			map = Basemap(projection=projection,lat_0=lat_0,lon_0 = lon_0,resolution = 'l',area_thresh=1000.0,celestial=True,ax = ax)
			index_,centor = index_list[n],centor_list[n]
			
			if source:
				ra, dec = map(source.ra.value, source.dec.value)
				map.plot(ra, dec, '*', color='#f36c21', markersize=20.)
			if points:
				ra, dec = map(points.ra.value, points.dec.value)
				map.plot(ra, dec, '*', color='#c7a252', markersize=20.)
			if show_bodies and sc_pos is not None:
				if projection in ['moll']:
					
					postion, r, lon, lat = earth_points_list[n]
					lon_lis, lat_lis = get_poly(postion, r, lon, lat, pole, lon_0)
					for i in range(len(lon_lis)):
						x, y = map(lon_lis[i], lat_lis[i])
						earth = Polygon(list(zip(x, y)), facecolor='#90d7ec',
						                edgecolor='#90d7ec',
						                linewidth=0, alpha=1)
						ax.add_patch(earth)
				else:
					postion, r, lon, lat = earth_points_list[n]
					lon, lat = map(lon, lat)
					earth = Polygon(list(zip(lon, lat)), facecolor='#90d7ec', edgecolor='#90d7ec',
					                linewidth=2, alpha=1)
					ax.add_patch(earth)
				if time is not None:
					earth_r = get_body_barycentric('earth',time[n])
					moon_r = get_body_barycentric('moon',time[n])
					r_e_m = moon_r - earth_r
					r = sc_pos[n] - np.array([r_e_m.x.value,r_e_m.y.value,r_e_m.z.value])*u.km
					moon_point_d = cartesian_to_spherical(-r[0],-r[1],-r[2])
					moon_ra,moon_dec = moon_point_d[2].deg,moon_point_d[1].deg
					moon_point = SkyCoord(moon_ra,moon_dec,frame='icrs', unit='deg')
					moon_ra,moon_dec = map(moon_point.ra.deg,moon_point.dec.deg)
					map.plot(moon_ra,moon_dec,'o',color = '#72777b',markersize = 20)
					plt.text(moon_ra,moon_dec-800000,'moon',size = 20)
				if show_bodies and time is not None:
					tmp_sun = get_sun(time[n])
					sun_position = SkyCoord(tmp_sun.ra.deg,tmp_sun.dec.deg,unit='deg', frame='icrs')
					sun_ra,sun_dec = map(sun_position.ra.value,sun_position.dec.value)
					map.plot(sun_ra,sun_dec ,'o',color = '#ffd400',markersize = 40)
					plt.text(sun_ra-550000,sun_dec-200000,'sun',size = 20)
			fovs = get_fov(centor,radius)
			if projection in ['moll']:
				for i,v in enumerate(index_):
					r,ra,dec = fovs[i]
					lon_lis, lat_lis = get_poly(centor[i],r,ra,dec,pole,lon_0)
					for ij in range(len(lon_lis)):
						x,y = map(lon_lis[ij],lat_lis[ij])
						detec = Polygon(list(zip(x,y)),facecolor=color_list[v],edgecolor=color_list[v],linewidth=2, alpha=0.5)
						ax.add_patch(detec)
					ra_x,dec_y = map(centor[i].ra.value+2.5,centor[i].dec.value-1)
					plt.text(ra_x, dec_y,str(name_list[v]), color=color_list[v], size=22)
			else:
				for i,v in enumerate(index_):
					r,ra,dec = fovs[i]
					detec = Polygon(list(zip(ra,dec)),facecolor=color_list[v],edgecolor=color_list[v],linewidth=2, alpha=0.5)
					ax.add_patch(detec)
					ra_x,dec_y = map(centor[i].ra.value+2.5,centor[i].dec.value-1)
					plt.text(ra_x, dec_y,str(name_list[v]), color=color_list[v], size=22)
			
			
			
			if projection == 'moll':
				az1 = np.arange(0, 360, 30)
				zen1 = np.zeros(az1.size) + 2
				azname = []
				for i in az1:
					azname.append(r'${\/%s\/^{\circ}}$' % str(i))
				x1, y1 = map(az1, zen1)
				for index1, value in enumerate(az1):
					plt.text(x1[index1], y1[index1], azname[index1], size=20)
			map.drawmeridians(np.arange(0, 360, 30),dashes=[1,0],color='#d9d6c3')
			map.drawparallels(np.arange(-90, 90, 15),dashes=[1,0], labels=[1,0,0,1], color='#d9d6c3',size = 20)
			map.drawmapboundary(fill_color='#f6f5ec')
			n = n + 1
			return mplfig_to_npimage(fig)
		animation = VideoClip(make_frame, duration=duration)
		#animation.write_videofile(savevdir, fps=1/dt,codec = 'h264')
		animation.write_videofile(savevdir, fps=1/dt)
		animation.close()
		

def get_fov(conter,radius = 10.):
		fov_point_list = []
		for conter_i in conter:
			poly = SphericalPolygon.from_cone(conter_i.ra.value,conter_i.dec.value,radius,steps=100)
			x,y = [p for p in poly.to_radec()][0]
			fov_point_list.append([radius,x,y])
		return fov_point_list

def loop_data(rad):
	rad = np.array(rad)
	rad[rad < 0] = rad[rad < 0] + 360
	rad[rad >= 360] = rad[rad >= 360] - 360
	return rad

def tranlation_lon_0(rad, lon_0):
	rad = np.array(rad)
	rad = rad - lon_0
	rad = loop_data(rad)
	return rad

def judge(x,n,xxdd,xxdd180):
	d_ra = np.abs(get_ra(x,n))
	case1 = (x>xxdd)&(xxdd>180)&(n<=xxdd)&(n>xxdd180)&(d_ra<90)
	case2 = (n>xxdd)&(xxdd>180)&(x<=xxdd)&(x>xxdd180)&(d_ra<90)
	case3 = (x<xxdd180)&(x>=xxdd)&(xxdd<=180)&(n<xxdd)&(d_ra<90)
	case4 = (n<xxdd180)&(n>=xxdd)&(xxdd<=180)&(x<xxdd)&(d_ra<90)
	return case1|case2|case3|case4

def get_ra(x,x0):
	xxx = np.arange(-5,6,1)
	xar = xxx*360+x
	dx = xar - x0
	dxabs = np.abs(dx)
	return dx[np.argmin(dxabs)]

def get_poly(position, radius, x, y, pole, lon_0):
	x = np.array(x)
	y = np.array(y)
	center_ra = position.ra.value
	xx_dd = tranlation_lon_0([lon_0 - 180], center_ra - 180)[0]
	x = tranlation_lon_0(x, center_ra - 180)
	inde = pole.separation(position) <= radius * u.degree
	pole_in = pole[inde]
	add_x = np.zeros(100) + xx_dd
	add_y = np.linspace(-89.9999, 89.9999, 100)
	if pole_in.size > 0:
		xx_dd180 = loop_data([xx_dd+180])[0]
		#print('pole_in.size > 0')
		if pole_in.dec.deg == 90.:
			#print('pole_in.dec.deg == 90.')
			x_n = []
			y_n = []
			n0 = x[0]
			y0 = y[0]
			for i in range(1, len(x)):
				
				if judge(x[i],n0,xx_dd,xx_dd180) == False :
					x_n.append(x[i])
					y_n.append(y[i])
					n0 = x[i]
					y0 = y[i]
				else:
					#print(x[i],n0,xx_dd,xx_dd180)
					indexs = np.where(add_y >= y0)[0]
					add_y_a = add_y[indexs]
					add_x_a = add_x[indexs]
					indexssort = np.argsort(add_y_a)
					ds = get_ra(n0,x[i])
					n_v = ds / np.abs(ds)
					add_y_a = add_y_a[indexssort]
					add_x_a = add_x_a[indexssort] + 0.1*n_v
					for i in range(len(add_y_a)):
						x_n.append(add_x_a[i])
						y_n.append(add_y_a[i])
					indexssort = np.argsort(-add_y_a)
					add_y_a = add_y_a[indexssort]
					
					add_x_a = add_x_a[indexssort] - 0.1*n_v
					for i in range(len(add_y_a)):
						x_n.append(add_x_a[i])
						y_n.append(add_y_a[i])
					n0 = x[i]
					y0 = y[i]
			x = [tranlation_lon_0(x_n, -center_ra + 180)]
			y = [np.array(y_n)]
		else:
			#print('pole_in.dec.deg == 90. else')
			x_n = []
			y_n = []
			n0 = x[0]
			y0 = y[0]
			
			for i in range(1, len(x)):
				if judge(x[i],n0,xx_dd,xx_dd180) == False:
					
					x_n.append(x[i])
					y_n.append(y[i])
					n0 = x[i]
					y0 = y[i]
				else:
					#print('dddddd',x[i],n0,xx_dd,xx_dd180)
					#print(i)
					indexs = np.where(add_y <= y0)[0]
					add_y_a = add_y[indexs]
					#add_x_a = np.zeros(add_y_a.size) + n0
					add_x_a = add_x[indexs]
					indexssort = np.argsort(-add_y_a)
					add_y_a = add_y_a[indexssort]
					ds = get_ra(n0,x[i])
					n_v = ds / np.abs(ds)
					add_x_a = add_x_a[indexssort]+0.1*n_v
					
					#print(n0,x[i],xx_dd)
					#print(n0 - xx_dd)
					for i in range(len(add_y_a)):
						#pass
						x_n.append(add_x_a[i])
						y_n.append(add_y_a[i])
					indexssort = np.argsort(add_y_a)
					add_y_a = add_y_a[indexssort]
					
					add_x_a = add_x_a[indexssort] - 0.1*n_v
				
					for i in range(len(add_y_a)):
						
						x_n.append(add_x_a[i])
						y_n.append(add_y_a[i])
					n0 = x[i]
					y0 = y[i]

			x = [tranlation_lon_0(x_n, -center_ra + 180)]
			y = [np.array(y_n)]
	else:
		#print('pole_in.size > 0 else')
		x_1 = []
		y_1 = []
		x_2 = []
		y_2 = []
		
		if len(x[x > xx_dd]) > len(x[x <= xx_dd]):
			#print('len(x[x > xx_dd]) > len(x[x <= xx_dd])')
			if len(x[x <= xx_dd]) > 0:
				#print('len(x[x <= xx_dd]) > 0')
				# x大于xx_dd
				for i in range(len(x)):
					# x小于xx_dd一律等于xx_dd+0.1
					if x[i] <= xx_dd:
						x_1.append(xx_dd + 0.1)
						y_1.append(y[i])
					else:
						x_1.append(x[i])
						y_1.append(y[i])
				
				aa = y[x <= xx_dd]
				amax = aa.max()
				amin = aa.min()
				for i in range(len(x)):
					# x大于xx_dd一律等于xx_dd-0.1
					if x[i] >= xx_dd:
						x_2.append(xx_dd - 0.1)
						if y[i] > amax:
							y_2.append(amax)
						elif y[i] < amin:
							y_2.append(amin)
						else:
							y_2.append(y[i])
					else:
						x_2.append(x[i])
						if y[i] > amax:
							y_2.append(amax)
						elif y[i] < amin:
							y_2.append(amin)
						else:
							y_2.append(y[i])
				x = [tranlation_lon_0(x_1, -center_ra + 180),
				     tranlation_lon_0(x_2, -center_ra + 180)]
				y = [np.array(y_1), np.array(y_2)]
			else:
				x = [tranlation_lon_0(x, -center_ra + 180)]
				y = [y]
		else:
			#print('len(x[x > xx_dd]) > len(x[x <= xx_dd]) else')
			if len(x[x > xx_dd]) > 0:
				#print('len(x[x > xx_dd]) > 0')
				# 小于xx_dd为主
				for i in range(len(x)):
					if x[i] >= xx_dd:
						x_1.append(xx_dd - 0.1)
						y_1.append(y[i])
					else:
						x_1.append(x[i])
						y_1.append(y[i])
				
				aa = y[x >= xx_dd]
				amax = aa.max()
				amin = aa.min()
				for i in range(len(x)):
					if x[i] <= xx_dd:
						x_2.append(xx_dd + 0.1)
						if y[i] > amax:
							y_2.append(amax)
						elif y[i] < amin:
							y_2.append(amin)
						else:
							y_2.append(y[i])
					else:
						x_2.append(x[i])
						if y[i] > amax:
							y_2.append(amax)
						elif y[i] < amin:
							y_2.append(amin)
						else:
							y_2.append(y[i])
				x = [tranlation_lon_0(x_1, -center_ra + 180),
				     tranlation_lon_0(x_2, -center_ra + 180)]
				y = [np.array(y_1), np.array(y_2)]
			else:
				#print('len(x[x > xx_dd]) > 0 else')
				x = [tranlation_lon_0(x, -center_ra + 180)]
				y = [y]
	return x, y
