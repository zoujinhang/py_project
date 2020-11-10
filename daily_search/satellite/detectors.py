

import numpy as np
from astropy.coordinates import spherical_to_cartesian

class Detectors(object):

	def __init__(self,local_az=None,local_zen=None,namelist = None,eff_angle=None):

		if local_az is None or local_zen is None or namelist is None or eff_angle is None:

			self.name = np.array(['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb','b0','b1'])
			self.detetor_vector = np.array(
				[[0.2446677589, 0.2523893824, 0.9361823057],
				 [0.5017318971, 0.5036621127, 0.7032706462],
				 [0.5233876659, 0.8520868147, -0.0036651682],
				 [0.5009495177, -0.5032279093, 0.7041386753],
				 [0.5468267487, -0.8372325378, -0.0047123847],
				 [0.9982910766, 0.0584352143, 0.0005236008],
				 [-0.2471260191, -0.2465229020, 0.9370993606],
				 [-0.5135631636, -0.5067957667, 0.6923950822],
				 [-0.5503349679, -0.8349438131, 0.0005235846],
				 [-0.5064476761, 0.5030998708, 0.7002865795],
				 [-0.5552650628, 0.8316411478, -0.0073303046],
				 [-0.9978547710, -0.0652279514, -0.0055850266],
				 [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
			)
			self.eff_angle = np.array([60,60,60,60,60,60,60,60,60,60,60,60,60,90,90])

		else:
			self.name = np.array(namelist)
			self.eff_angle = np.array(eff_angle)
			self.detetor_vector = np.array(spherical_to_cartesian(1,local_az,local_zen)).T



	def __call__(self,name):

		n_index = np.where(self.name == name)[0]
		detector_vector = self.detetor_vector[n_index]

		return detector_vector[0]

	def get_eff_angle(self,name):

		n_index = np.where(self.name == name)[0]
		eff_angle = self.eff_angle[n_index]
		return eff_angle[0]

	def get_names(self):
		return self.name

