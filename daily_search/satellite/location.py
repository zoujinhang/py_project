
from astropy.coordinates import cartesian_to_spherical,SkyCoord
import astropy.units as u
import numpy as np
from ..file import findfile,readcol
import sys
import os
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt



file_name_list = ['alocdat_comp.dat',
		  'earth_points.dat',
		  'locrates_1deg_50_300_hard_n.dat',
		  'locrates_1deg_50_300_normal_n.dat',
		  'locrates_1deg_50_300_soft_n.dat',
		  'locrates_1deg_5_50_soft.dat']


data_name_list = ['scatterdata','earthpoints','h','n','s','sl']


paths = sys.path
location_data = {}

for i in range(len(data_name_list)):

	for path in paths:
		fand = path+'/daily_search/satellite/responce_data/'
		if os.path.exists(fand):
			sonamelist = findfile(fand,file_name_list[i])
			if len(sonamelist)>0:

				location_data[data_name_list[i]] = fand+sonamelist[0]
				print('the C lib link is ',fand+sonamelist[0])
				break
			else:
				print('we do not find the file ',file_name_list[i])
				raise FileNotFoundError

earthpoints = np.array(readcol(location_data['earthpoints']))
scatterdata = np.array(readcol(location_data['scatterdata'])).T.reshape((19,236,16))
scatterdata = scatterdata[:,:,::-1]
'''
scatterdata = scatterdata[:,:,::-1]

fig = plt.figure(figsize = (20,50))
ax = fig.gca(projection = '3d')
for i in range(16):
	for j in range(19):
		for k in range(236):
			ax.scatter(i,j,k,marker = ',',c = (scatterdata[i,j,k]-scatterdata.min())/(scatterdata.max()-scatterdata.min()),
				   alpha=(scatterdata[i,j,k]-scatterdata.min())/(scatterdata.max()-scatterdata.min()))
fig.savefig('/home/laojin/my_lat/location/A_data.png')
plt.close(fig)
'''

case_n = {
	#'case1':(np.array(readcol(location_data['sl'])),np.array([-2.0, -3.4, 70., 10.])),
	'case1':(np.array(readcol(location_data['sl'])),None),
	'case2':(np.array(readcol(location_data['s'])),np.array([-2.0, -3.4, 70., 10.])),
	'case3':(np.array(readcol(location_data['n'])),np.array([-1.0, -2.3, 230., 10.])),
	'case4':(np.array(readcol(location_data['h'])),np.array([0., -1.5, 1000., 10.])),
	}

def crossprod(x,y):

	return np.vstack([x[1]*y[2] - x[2]*y[1],
			  x[2]*y[0] - x[0]*y[2],
			  x[0]*y[1] - x[1]*y[0]])


class Locate(object):

	def __init__(self,geometry):
		'''

		'''

		self.geometry = geometry
		self.earthpoints = earthpoints
		self.scatterdata = scatterdata
		self.tenergies = np.array([9.88152, 21.9039, 30.6248, 39.6809, 53.7016, 72.9625, 97.6607, 122.844, 163.847,
					   231.738, 316.693, 424.295, 587.606, 741.422, 1096.22, 1813.54, 2749.15])
		self.case_n = case_n
		self.detector = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
		self.grid_points = np.linspace(0,180,19)/180*np.pi

	def condidate(self,t,rate,bs_rate,case,cenergies):
		'''

		'''
		pos_ = self.geometry.get_pos(t)
		qsj = self.geometry.get_qsj(t)
		loc_pos = self.geometry.inv_transform_frame_one(pos_,qsj)
		loc_pos_mag = np.sqrt((loc_pos**2).sum())
		loc_pos = loc_pos/loc_pos_mag
		earth_r = self.geometry.get_earth_point(t)[-1]
		loctable_entries = self.get_entries(case,loc_pos,cenergies,earth_r)
		az = loctable_entries[0]/60/180*np.pi
		el = loctable_entries[1]/60/180*np.pi
		r_cart = np.vstack([np.sin(el)*np.cos(az),
				    np.sin(el)*np.sin(az),
				    np.cos(el)]).T

		entries = (loctable_entries[2:]).T
		fnorm = ((entries*(rate-bs_rate)/rate).sum(axis = 1))/(entries**2/rate).sum(axis = 1)
		ccc = (fnorm*(entries.T)).T
		chi2 = ((rate-bs_rate-ccc)**2/(bs_rate+ccc)).sum(axis = 1)
		chi2[chi2<0] = 99999
		sort_index = np.argsort(chi2)
		gindex = sort_index[0]
		gindex2 = sort_index[1]
		tran_cart = self.geometry.transform_frame_more(r_cart,qsj)
		position_cart = cartesian_to_spherical(tran_cart[:,0],tran_cart[:,1],tran_cart[:,2])

		good_cart = r_cart[gindex]
		good_cart = self.geometry.transform_frame_one(good_cart,qsj)
		position = cartesian_to_spherical(good_cart[0],good_cart[1],good_cart[2])
		error_determined = False
		loc_reliable = True
		loc_err_vsmall = False
		offset_chi2_delta = 0.025
		if (chi2[gindex2]-chi2[gindex] > 2.3):
			loc_err_vsmall = True
			loc_err = 1.
		## test code for fiducial intensity check for chi2 reliability
		d_chi2 = chi2-chi2[gindex]
		print(np.sort(d_chi2))
		while error_determined == False and loc_reliable and loc_err_vsmall==False:

			index_chi = np.where((d_chi2>= 9.21-offset_chi2_delta)&(d_chi2<=9.21+offset_chi2_delta))[0]
			if len(index_chi)>1:
				error_determined = True
				xyz_position = SkyCoord(x=r_cart[index_chi,0],y=r_cart[index_chi,1],z=r_cart[index_chi,2],frame='icrs',representation='cartesian')
				good_position = SkyCoord(x=r_cart[gindex,0],y=r_cart[gindex,1],z=r_cart[gindex,2],frame='icrs',representation='cartesian')
				sep_ = xyz_position.separation(good_position).deg
				print('what',sep_)
				loc_err = np.mean(sep_)
			if (offset_chi2_delta >= 9.21):
				loc_err=50.
				loc_reliable = False
				error_determined = True
				print('location unreliable')
			offset_chi2_delta=offset_chi2_delta+0.05

		return position[2].deg,position[1].deg,loc_err,position_cart[2].deg,position_cart[1].deg,chi2




	def get_spec(self,param,tenergies,erange):


		mid_enegy = np.sqrt(tenergies[1:]*tenergies[:-1])
		erange_index = np.where((mid_enegy>=erange[0])&(mid_enegy<=erange[1]))

		d_energy = tenergies[1:] - tenergies[:-1]
		zer = np.zeros(len(d_energy))
		zer[erange_index]=1
		alpha = param[0]
		beta = param[1]
		epeak= param[2]
		fnorm = param[3]
		if(alpha < -1.9):
			alpha=-1.9
		e0=epeak/(2.+alpha)
		eb=(alpha-beta)*e0
		spec = np.zeros(len(mid_enegy))
		index_1 = np.where(mid_enegy<=eb)[0]
		spec[index_1] = (mid_enegy[index_1]/100)**alpha * np.exp(-mid_enegy[index_1]/e0)
		index_2 = np.where(mid_enegy>eb)[0]
		spec[index_2] = ((alpha-beta)*e0/100.)**(alpha-beta)* np.exp(beta-alpha)*(mid_enegy[i]/100.)**beta
		spec = spec * d_energy
		f=(spec[erange_index]).sum()
		spec = spec*zer
		return spec*fnorm/f


	def get_entries(self,case,loc_pos,cenergies,earth_r):


		geom_fac_front=126*50/(2025.*0.8)
		geom_fac_back=126*50/(2025.*0.8)

		entries , spec = self.case_n[case]


		az = entries[0]/60/180*np.pi    #rad
		el = entries[1]/60/180*np.pi    #rad
		r_cart = np.vstack([np.sin(el)*np.cos(az),
				    np.sin(el)*np.sin(az),
				    np.cos(el)]).T

		xyz_position = SkyCoord(x=r_cart[:,0],y=r_cart[:,1],z=r_cart[:,2],frame='icrs',representation='cartesian')
		earthp = SkyCoord(x=loc_pos[0],y=loc_pos[1],z=loc_pos[2],frame='icrs',representation='cartesian')
		seq = xyz_position.separation(earthp) > earth_r*u.degree
		entries = entries[:,seq]
		r_cart = r_cart[seq]
		if spec is None:
		#if True:
			return entries
		scat_spec = self.get_spec(spec, self.tenergies, cenergies)
		atm_scattered_rates = np.zeros((len(self.detector), len(r_cart)))
		#eearthpoints = np.zeros((3,236))
		#heearthpoints = np.zeros((3,236))
		for direction in [1,-1]:
			geom_fac = geom_fac_front
			scat_geom_fac = 1
			if direction != 1:
				#geom_fac = geom_fac_back
				#scat_geom_fac = 0
				continue

			for dete_index,dete in enumerate(self.detector):
				detec_ver = direction*self.geometry.detectors(dete)
				detec_ver_mag = np.sqrt((detec_ver ** 2).sum())
				detec_ver = detec_ver / detec_ver_mag
				cos_sita = (detec_ver * loc_pos).sum()
				az_ = np.zeros(len(r_cart))
				if cos_sita > 1:
					cos_sita = 1
				if cos_sita < -1:
					cos_sita = -1
				angle_rad = np.arccos(cos_sita)
				elev_idet = np.pi * 0.5 - angle_rad
				if elev_idet < -1.48353:
					elev_idet = -1.48353
				if elev_idet > 1.48353:
					elev_idet = 1.48353
				gperp_idet = -loc_pos - cos_sita * detec_ver
				gperp_idet_mag = np.sqrt((gperp_idet**2).sum())

				c = crossprod(r_cart.T,detec_ver)
				a = (crossprod(c,detec_ver)).T
				amag = np.sqrt((a*a).sum(axis = 1))
				if gperp_idet_mag>1e-5:
					amag_index = np.where(amag > 1e-5)[0]
					if len(amag_index)>0:
						x = (gperp_idet*a[amag_index]).sum(axis=1)/(gperp_idet_mag*amag[amag_index])
						x[x>1]=1
						x[x<-1] = -1
						az_[amag_index] = np.arccos(x)
				eindex = 0
				while self.earthpoints[2, eindex] > elev_idet:
					eindex = eindex + 1
				heindex = eindex - 1
				bdangle = np.arccos((detec_ver * r_cart).sum(axis=1))
				#print('--------------------------------------')
				'''
				eearthpoints[:,:]=99.
				heearthpoints[:,:]=99.
				for i in range(236):
					if self.earthpoints[2,i] == self.earthpoints[2, eindex]:
						eearthpoints[:,i] = self.earthpoints[:, i]
					if self.earthpoints[2,i] == self.earthpoints[2, heindex]:
						heearthpoints[:,i] = self.earthpoints[:, i]
				'''
				earthpoints_index = self.earthpoints[2] == self.earthpoints[2, eindex]
				eearthpoints = self.earthpoints[:,earthpoints_index]

				earthpoints_index = self.earthpoints[2] == self.earthpoints[2, heindex]
				heearthpoints = self.earthpoints[:,earthpoints_index]

				#print(eearthpoints)
				#print(eearthpoints1)
				#print('--------------------------------------')
				scattered_rates = np.zeros(16)
				for r_c_index,r_cart_i in enumerate(r_cart):
					#print(r_c_index)
					'''
					c = crossprod(r_cart_i,detec_ver)
					a = crossprod(c,detec_ver)
					amag = np.sqrt((a*a).sum())
					if gperp_idet_mag>1e-5 and amag > 1e-5:
						x = (gperp_idet*a/(gperp_idet_mag*amag)).sum()
						if x > 1.:
							x=1.
						if x <-1. :
							x=-1
						az_ = np.arccos(x)
					else:
						az_ = 0
					'''
					aindex = np.argmin(np.abs(eearthpoints[1]-az_[r_c_index]))
					#print('aindex1',aindex)
					aindex = int(eearthpoints[0,aindex])-1
					#aindex2 = np.argmin(np.abs(eearthpoints1[1]-az_[r_c_index]))
					#print('aindex1',aindex)
					#aindex2 = int(eearthpoints1[0,aindex2])-1
					haindex = np.argmin(np.abs(heearthpoints[1]-az_[r_c_index]))
					haindex = int(heearthpoints[0,haindex])-1

					felev = 0
					if self.earthpoints[2,haindex]-self.earthpoints[2,aindex] > 0.001:
						felev=(elev_idet-self.earthpoints[2,aindex])/(self.earthpoints[2,haindex]-self.earthpoints[2,aindex])
					#bdanglei = np.arccos((detec_ver * r_cart_i).sum())  # !!!!

					scattered_rates[:] = 0.
					for i in range(16):
						#print('scat',self.scatterdata[:,aindex,i])
						r1 = np.interp(bdangle[r_c_index],self.grid_points,self.scatterdata[:,aindex,i])
						r2= np.interp(bdangle[r_c_index],self.grid_points,self.scatterdata[:,haindex,i])
						scattered_rates[i] = scat_geom_fac*((1.-felev)*r1+felev*r2)
					#print('scatered_rates',scattered_rates)
					atm_scattered_rates[dete_index,r_c_index] = (geom_fac*scattered_rates*scat_spec).sum()+atm_scattered_rates[dete_index,r_c_index]
		entries[2:] = entries[2:]+atm_scattered_rates
		return entries









