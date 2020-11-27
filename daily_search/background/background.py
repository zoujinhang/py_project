

from .baseline import TD_baseline
from scipy.interpolate import interp1d


def get_background_f(t,rate,**kwargs):
	#ta = np.zeros(len(t))
	bs = TD_baseline(t,rate,**kwargs)
	dt = t[1]-t[0]
	#print('dt_get_background',dt)
	ta = t + 0.
	ta[0] = ta[0]-dt
	ta[-1] = ta[-1] + dt
	bs_f = interp1d(ta,bs,kind = 'slinear') #

	return bs_f








