
from Fermi_tool.daily import Database,search_candidates,Sources,track,Plot_track,Save_track,Save_search,time_slic,Plot_serch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
from Data_analysis.geometry import Geometry
from multiprocessing import Pool


savedir = '/home/laojin/my_lat/daily_text/sgr1/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)


#topdir = '/media/laojin/Elements/daily/'
#topdir = '/media/laojin/TOSHIBA_EXT/daily/'
topdir = '/home/laojin/daily_data/'
timestart = '2020-04-28T00:00:00'
timestop = '2020-04-28T00:59:00'

#timestart = '2013-04-27T07:45:00'
#timestop = '2013-04-27T07:53:00'
ra = [293.729]
dec = [21.3864]
name = ['SGRJ1935']
#ra = [170.6667,299,344.08]
#dec = [48.0500,7,-60.0]
#name = ['bn130427324','bn200219413','bn200219317']
# ----------------------------------------------------
# Create the sources
mysources = Sources(positions=[ra, dec], names=name, range_=5)
# ---------------------------------------------------_
t_sl = time_slic(timestart,timestop,H=3)
print('The time interval, H =',3)
print(t_sl)
print(str(t_sl[0][0]))
t_savedir = savedir
print(t_savedir)
if os.path.exists(t_savedir) ==False:
	os.makedirs(t_savedir)
s_savedir = t_savedir + 'Z_no_direction_trig/'
if os.path.exists(s_savedir) == False:
	os.makedirs(s_savedir)

#for i,t_s in enumerate(t_sl):

def run(t_s):

	t_start, t_stop = t_s
	time_markker = str(t_start)[:13]

	# ----------------------------------------------------
	# Create the geometry
	geometry = Geometry()
	# ----------------------------------------------------
	# Create the Database
	fermi = Database(topdir)
	# ----------------------------------------------------
	pos = fermi.get_poshist_data(t_start,t_stop)
	geometry.input_pose(quaternion=pos)
	detector = fermi.detector
	data = fermi.get_detector_data(t_start,t_stop)
	serch_result = search_candidates(data,detector,geometry)
	result = track(serch_result,geometry,mysources)
	plt_s = Plot_serch(serch_result,detector,geometry)
	

	bayes_size = plt_s.get_bayesian_responses_size()
	if bayes_size >0:
		for i in range(bayes_size):
			plt_s.plot_bayesian_responses(i)
			plt.savefig(s_savedir + 'A_bayes_'+time_markker+'_'+str(i)+'.png')
			plt.close()
	thres_size = plt_s.get_threshold_responses_size()
	if thres_size>0:
		for i in range(thres_size):
			
			plt_s.plot_threshold_responses(i)
			plt.savefig(s_savedir + 'B_threshold_'+time_markker+'_'+ str(i) + '.png')
			plt.close()
	
	ss = Save_search(serch_result,geometry)
	ss.save_all_responses(s_savedir + 'Z_save_total_'+time_markker+'.csv')
	ss.save_bayesian_responses(s_savedir + 'Z_save_bayesian_trig_'+time_markker+'.csv')
	ss.save_threshold_responses(s_savedir + 'Z_save_threshold_trig_'+time_markker+'.csv')
	
	
	plt_t = Plot_track(result, detector, geometry,sources=mysources)
	st = Save_track(result,geometry)
	#tirg_data = result['lc']
	#print(data['n0']['events'])
	for in_sn,sn in enumerate(mysources.names):


		sn_savedir = t_savedir + sn + '/'
		if os.path.exists(sn_savedir) == False:
			os.makedirs(sn_savedir)
		good_savedir_bayes = sn_savedir + 'A_bayesian_trig/'
		if os.path.exists(good_savedir_bayes) == False:
			os.makedirs(good_savedir_bayes)
		good_savedir = sn_savedir + 'A_simple_trig/'
		if os.path.exists(good_savedir) == False:
			os.makedirs(good_savedir)

		st.save_all_responses(sn,sn_savedir + 'B_save_total_'+time_markker+'.csv')
		st.save_bayesian_responses(sn,good_savedir_bayes + 'B_save_bayesian_trig_'+time_markker+'.csv')
		st.save_threshold_responses(sn,good_savedir + 'B_save_simple_trig_'+time_markker+'.csv')
		#plt.figure(constrained_layout=True, figsize=(20, 5*len(detector)))
		fig,axs = plt.subplots(nrows=len(detector)*2,figsize=(20, 5*len(detector)),constrained_layout=True)
		source = mysources.positions[in_sn]
		range_i = mysources.range
		plt_t.plot_one_source(sn,source,range_i,axs)
		fig.savefig(sn_savedir + 'B_'+sn+'_analysis_'+time_markker+'.png')
		plt.close(fig)
		

		
		bayes_size = plt_t.get_bayesian_responses_size(sn)
		if bayes_size >0:
			for i in range(bayes_size):
				plt_t.plot_bayesian_responses(sn,i,sky_map = True)
				plt.savefig(good_savedir_bayes + 'A_bayes_'+time_markker+'_'+str(i)+'.png')
				plt.close()

		thres_size = plt_t.get_threshold_responses_size(sn)
		if thres_size>0:
			for i in range(thres_size):
				plt_t.plot_threshold_responses(sn,i,sky_map = True)
				plt.savefig(good_savedir + 'A_threshold_'+time_markker+'_' + str(i) + '.png')
				plt.close()
		
		#all_size = plt_t.get_all_responses_size(sn)
		#if all_size>0:
		#	for i in range(all_size):
		#
		#		plt_t.plot_all_responses(sn,i)
		#		plt.savefig(sn_savedir+'B_all_'+str(i)+'.png')
		#		plt.close()

pool = Pool(2)     # Parallel computing
pool.map(run,list(t_sl))
pool.close()
pool.join()

