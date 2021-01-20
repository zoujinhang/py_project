
from daily_search import Database
import astropy.units as u
from daily_search.satellite import Geometry
from astropy.io import fits
import numpy as np
import os
from gbm_drm_gen.drmgen import DRMGen

datatop = '/media/laojin/My_Passport/Fermi_GBM_Dailiy_Data/'
savetop = '/home/laojin/my_lat/make_datafits/'
if os.path.exists(savetop) ==False:
	os.makedirs(savetop)

ra = 293.729
dec = 21.3864

outfile = savetop + 'A_lab_.fit'
outfile2 = savetop + 'A_lab_2.fit'
t_start = '2020-04-27T18:00:00'
t_stop = '2020-04-27T19:00:00'

det_name_lookup = {
    "n0": 0,
    "n1": 1,
    "n2": 2,
    "n3": 3,
    "n4": 4,
    "n5": 5,
    "n6": 6,
    "n7": 7,
    "n8": 8,
    "n9": 9,
    "na": 10,
    "nb": 11,
    "b0": 12,
    "b1": 13,
}

class DRMGentte(DRMGen):

	def __init__(self,quaternions,sc_pos,out_edge,deter,
		     mat_type=0,
		     occult=False,
		     ):
		self._det_number = det_name_lookup[deter]
		self._out_edge = out_edge
		self._matrix_type = mat_type
		self._occult = occult
		if self._det_number > 11:
			self._in_edge = np.array(
				[
				    100.000,
				    105.579,
				    111.470,
				    117.689,
				    124.255,
				    131.188,
				    138.507,
				    146.235,
				    154.394,
				    163.008,
				    172.103,
				    181.705,
				    191.843,
				    202.546,
				    213.847,
				    225.778,
				    238.375,
				    251.675,
				    265.716,
				    280.541,
				    296.194,
				    312.719,
				    330.167,
				    348.588,
				    368.036,
				    388.570,
				    410.250,
				    433.139,
				    457.305,
				    482.820,
				    509.757,
				    538.198,
				    568.226,
				    599.929,
				    633.401,
				    668.740,
				    706.052,
				    745.444,
				    787.035,
				    830.946,
				    877.307,
				    926.255,
				    977.933,
				    1032.49,
				    1090.10,
				    1150.92,
				    1215.13,
				    1282.93,
				    1354.51,
				    1430.08,
				    1509.87,
				    1594.11,
				    1683.05,
				    1776.95,
				    1876.09,
				    1980.77,
				    2091.28,
				    2207.96,
				    2331.15,
				    2461.21,
				    2598.53,
				    2743.51,
				    2896.58,
				    3058.18,
				    3228.81,
				    3408.95,
				    3599.15,
				    3799.96,
				    4011.97,
				    4235.81,
				    4472.14,
				    4721.65,
				    4985.09,
				    5263.22,
				    5556.87,
				    5866.90,
				    6194.24,
				    6539.83,
				    6904.71,
				    7289.95,
				    7696.67,
				    8126.09,
				    8579.47,
				    9058.15,
				    9563.53,
				    10097.1,
				    10660.5,
				    11255.2,
				    11883.2,
				    12546.2,
				    13246.2,
				    13985.2,
				    14765.5,
				    15589.3,
				    16459.1,
				    17377.4,
				    18346.9,
				    19370.6,
				    20451.3,
				    21592.4,
				    22797.1,
				    24069.0,
				    25411.8,
				    26829.7,
				    28326.6,
				    29907.0,
				    31575.6,
				    33337.3,
				    35197.3,
				    37161.0,
				    39234.4,
				    41423.4,
				    43734.5,
				    46174.6,
				    48750.8,
				    51470.7,
				    54342.5,
				    57374.4,
				    60575.5,
				    63955.2,
				    67523.4,
				    71290.7,
				    75268.2,
				    79467.7,
				    83901.5,
				    88582.6,
				    93524.9,
				    98742.9,
				    104252.0,
				    110069.0,
				    116210.0,
				    122693.0,
				    129539.0,
				    136766.0,
				    144397.0,
				    152453.0,
				    160959.0,
				    169939.0,
				    179421.0,
				    189431.0,
				    200000.0,
				],
				dtype=np.float32,)
		else:
			self._in_edge = np.array(
				[
				    5.00000,
				    5.34000,
				    5.70312,
				    6.09094,
				    6.50513,
				    6.94748,
				    7.41991,
				    7.92447,
				    8.46333,
				    9.03884,
				    9.65349,
				    10.3099,
				    11.0110,
				    11.7598,
				    12.5594,
				    13.4135,
				    14.3256,
				    15.2997,
				    16.3401,
				    17.4513,
				    18.6380,
				    19.9054,
				    21.2589,
				    22.7045,
				    24.2485,
				    25.8974,
				    27.6584,
				    29.5392,
				    31.5479,
				    33.6931,
				    35.9843,
				    38.4312,
				    41.0446,
				    43.8356,
				    46.8164,
				    50.0000,
				    53.4000,
				    57.0312,
				    60.9094,
				    65.0513,
				    69.4748,
				    74.1991,
				    79.2446,
				    84.6333,
				    90.3884,
				    96.5349,
				    103.099,
				    110.110,
				    117.598,
				    125.594,
				    134.135,
				    143.256,
				    152.997,
				    163.401,
				    174.513,
				    186.380,
				    199.054,
				    212.589,
				    227.045,
				    242.485,
				    258.974,
				    276.584,
				    295.392,
				    315.479,
				    336.931,
				    359.843,
				    384.312,
				    410.446,
				    438.356,
				    468.164,
				    500.000,
				    534.000,
				    570.312,
				    609.094,
				    650.512,
				    694.748,
				    741.991,
				    792.446,
				    846.333,
				    903.884,
				    965.349,
				    1030.99,
				    1101.10,
				    1175.98,
				    1255.94,
				    1341.35,
				    1432.56,
				    1529.97,
				    1634.01,
				    1745.13,
				    1863.80,
				    1990.54,
				    2125.89,
				    2270.45,
				    2424.85,
				    2589.74,
				    2765.84,
				    2953.92,
				    3154.79,
				    3369.31,
				    3598.43,
				    3843.12,
				    4104.46,
				    4383.56,
				    4681.65,
				    5000.00,
				    5340.00,
				    5703.12,
				    6090.94,
				    6505.12,
				    6947.48,
				    7419.91,
				    7924.46,
				    8463.33,
				    9038.84,
				    9653.49,
				    10309.9,
				    11011.0,
				    11759.8,
				    12559.4,
				    13413.5,
				    14325.6,
				    15299.7,
				    16340.1,
				    17451.3,
				    18637.9,
				    19905.3,
				    21258.9,
				    22704.5,
				    24248.5,
				    25897.3,
				    27658.4,
				    29539.2,
				    31547.8,
				    33693.1,
				    35984.3,
				    38431.2,
				    41044.6,
				    43835.6,
				    46816.4,
				    50000.0,
				],
				dtype=np.float32,)
		super(DRMGentte, self).__init__(quaternions=quaternions,
						sc_pos=sc_pos,
						det_number=self._det_number,
						ebin_edge_in=self._in_edge,
						mat_type=self._matrix_type,
						ebin_edge_out=self._out_edge,
						occult=self._occult,)

fermi_data = Database(datatop)

data = fermi_data.get_detector_data(t_start,t_stop)
pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)

geometry = Geometry(pd_pos_data)
clock = geometry.clock
time_met = clock.utc_to_met('2020-04-27T18:33:05.850')
qsj = geometry.get_qsj(time_met)
pos = -geometry.get_pos(time_met)*u.m

for detei in data.keys():
	ni = data[detei]
	ch_E = ni['ch_E']
	out_edge = np.zeros(129, dtype=np.float32)
	out_edge[:-1] = ch_E["E_MIN"].values
	out_edge[-1] = ch_E["E_MAX"].values[-1]
	detei_i = DRMGentte(qsj,pos,out_edge,detei)
	savename = savetop + 'B_'
	detei_i.to_fits(ra,dec,savename,overwrite=True)

'''
name = data.keys()
ni = data['n0']
ch_E = ni['ch_E']
t = ni['events']['TIME'].values
ch = ni['events']['PHA'].values

pd_pos_data = fermi_data.get_poshist_data(t_start,t_stop)

header0=fits.Header()
header0.append(('creator', 'Zou', 'The name who created this PHA file'))
header0.append(('telescop', 'Fermi', 'Name of mission/satellite'))
header0.append(('TSTART', 0,'Observation start time'))
header0.append(('TSTOP', 0,'Observation stop time'))
header0.append(('TRIGTIME',0,'Trigger time relative to MJDREF, double precisi'))
header0.append(('TRIGA', 0,'Trigger start'))
header0.append(('TRIGB', 0,'Trigger stop'))

hdu0=fits.PrimaryHDU(header=header0)
a1 = np.arange(128)
col1 = fits.Column(name='CHANNEL', format='1I', array=ch_E['CHANNEL'].values)                              #创建列
col2 = fits.Column(name='E_MIN', format='1E', unit='keV', array=ch_E['E_MIN'].values)
col3 = fits.Column(name='E_MAX', format='1E', unit='keV', array=ch_E['E_MAX'].values)
hdu1 = fits.BinTableHDU.from_columns([col1, col2, col3])

col4 = fits.Column(name='TIME', format='1D', unit='s', array=t)                              #创建列
col5 = fits.Column(name='PHA', format='1I', unit='none', array=ch)
hdu2 = fits.BinTableHDU.from_columns([col4, col5])

hdul = fits.HDUList([hdu0, hdu1,hdu2])

if(os.path.exists(outfile)):
	os.remove(outfile)  # 删除旧版本

hdul.writeto(outfile)
hdul.close()

header0=fits.Header()
header0.append(('creator', 'Zou', 'The name who created this PHA file'))
header0.append(('telescop', 'Fermi', 'Name of mission/satellite'))
header0.append(('TSTART', 0,'Observation start time'))
header0.append(('TSTOP', 0,'Observation stop time'))
header0.append(('TRIGTIME',0,'Trigger time relative to MJDREF, double precisi'))
header0.append(('TRIGA', 0,'Trigger start'))
header0.append(('TRIGB', 0,'Trigger stop'))

hdu0=fits.PrimaryHDU(header=header0)
a1 = np.arange(128)
col_name = ['SCLK_UTC','QSJ_1','QSJ_2','QSJ_3','QSJ_4','POS_X','POS_Y','POS_Z','SC_LAT','SC_LON']
unit = ['s',None,None,None,None,'m','m','m','deg','deg']
formatlist = ['1D','1D','1D','1D','1D','1E','1E','1E','1E','1E']
collist = []
for i in range(len(col_name)):
	coli = fits.Column(name=col_name[i], format=formatlist[i], unit=unit[i],array=pd_pos_data[col_name[i]].values)
	collist.append(coli)
hdu1 = fits.BinTableHDU.from_columns(collist)

hdul = fits.HDUList([hdu0, hdu1])

if(os.path.exists(outfile2)):
	os.remove(outfile2)  # 删除旧版本

hdul.writeto(outfile2)
hdul.close()

'''