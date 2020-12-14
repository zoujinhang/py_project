
import numpy as np
from astropy.coordinates import spherical_to_cartesian
import astropy.units as u


class Detectors(object):

	def __init__(self,name='a'):
		az = np.array([0,40,40,40,40,40,40,73.5,73.5,73.5,73.5,73.5,73.5,73.5,73.5,73.5,73.5,90,90,90,90,90,90,90,90])*u.deg
		zen = np.array([0,332.17,278.67,225.17,152.17,98.67,45.17,350.5,314.5,278.5,242.5,206.5,170.5,134.5,98.5,62.5,26.5,0,305,270,215,180,125,90,35])*u.deg
		self.detetor_vector = np.array(spherical_to_cartesian(1,90*u.deg-az,zen)).T
		name_list = []
		for i in range(len(az)):
			name_list.append(name+str(i))
		self.name = np.array(name_list)
		self.eff_angle = np.array([60]*25)

	def __call__(self, name):

		n_index = np.where(self.name == name)[0]
		detector_vector = self.detetor_vector[n_index]
		return detector_vector[0]

	def get_eff_angle(self,name):

		n_index = np.where(self.name == name)[0]
		eff_angle = self.eff_angle[n_index]
		return eff_angle[0]

	def get_names(self):
		return self.name