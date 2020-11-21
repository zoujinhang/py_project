
from .clock import Clock
from .background import get_background_f
from .satellite import Detectors,Geometry,Locate
from .perception import Event


def search_trig(data):

	name = data.keys()
	for deteri in name:

		ni = data[deteri]['events']
		ch_E = data[deteri]['ch_E']







