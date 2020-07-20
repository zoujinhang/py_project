
from Fermi_tool.daily import Database,search_candidates,Sources,track,Plot_track,Save_track,Save_search,time_slic,Plot_serch
import numpy as np
import matplotlib.pyplot as plt
import os
from Data_analysis.geometry import Geometry
from multiprocessing import Pool


savedir = '/home/laojin/my_lat/daily_text/0/'
if os.path.exists(savedir) ==False:
	os.makedirs(savedir)


topdir = '/media/laojin/Elements/daily/'
timestart = '2013-04-27T00:00:00'
timestop = '2013-04-27T23:59:00'

#timestart = '2013-04-27T07:45:00'
#timestop = '2013-04-27T07:53:00'
ra = [170.6667,299,344.08]
dec = [48.0500,7,-60.0]
name = ['bn130427324','bn200219413','bn200219317']

# ----------------------------------------------------
# Create the sources
mysources = Sources(positions=[ra, dec], names=name, range_=10)
# ---------------------------------------------------_
t_sl = time_slic(timestart,timestop,H=3)
print('The time interval, H =',3)
print(t_sl)
print(str(t_sl[0][0]))


#for i,t_s in enumerate(t_sl):

def run(t_s):
	t_start, t_stop = t_s
	t_savedir = savedir +str(t_start)[:13] +'/'
	print(t_savedir)
	if os.path.exists(t_savedir) ==False:
		os.makedirs(t_savedir)
	
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
	
	s_savedir = t_savedir + 'good/'
	if os.path.exists(s_savedir) ==False:
		os.makedirs(s_savedir)
	bayes_size = plt_s.get_bayesian_responses_size()
	if bayes_size >0:
		for i in range(bayes_size):
			
			plt_s.plot_bayesian_responses(i)
			plt.savefig(s_savedir + 'A_bayes_'+str(i)+'.png')
			plt.close()
	thres_size = plt_s.get_threshold_responses_size()
	if thres_size>0:
		for i in range(thres_size):
			
			plt_s.plot_threshold_responses(i)
			plt.savefig(s_savedir + 'B_threshold_' + str(i) + '.png')
			plt.close()
	
	ss = Save_search(serch_result,geometry)
	ss.save_all_responses(t_savedir + 'Z_save_all.cvs')
	ss.save_bayesian_responses(t_savedir + 'Z_save_bayesian.cvs')
	ss.save_threshold_responses(t_savedir + 'Z_save_threshold.cvs')
	
	
	plt_t = Plot_track(result, detector, geometry)
	st = Save_track(result,geometry)
	#tirg_data = result['lc']
	#print(data['n0']['events'])
	for in_sn,sn in enumerate(mysources.names):
		sn_savedir = t_savedir + sn + '/'
		if os.path.exists(sn_savedir) == False:
			os.makedirs(sn_savedir)
		st.save_all_responses(sn,sn_savedir + 'A_save_all.cvs')
		st.save_bayesian_responses(sn,sn_savedir + 'A_save_bayesian.cvs')
		st.save_threshold_responses(sn,sn_savedir + 'A_save_threshold.cvs')
		#plt.figure(constrained_layout=True, figsize=(20, 5*len(detector)))
		fig,axs = plt.subplots(nrows=len(detector)*2,figsize=(20, 5*len(detector)),constrained_layout=True)
		source = mysources.positions[in_sn]
		range_i = mysources.range
		plt_t.plot_one_source(sn,source,range_i,axs)
		fig.savefig(t_savedir + 'C_'+sn+'_analysis.png')
		plt.close(fig)
		
		good_savedir = sn_savedir + 'good/'
		if os.path.exists(good_savedir) == False:
			os.makedirs(good_savedir)
		
		bayes_size = plt_t.get_bayesian_responses_size(sn)
		if bayes_size >0:
			for i in range(bayes_size):
				plt_t.plot_bayesian_responses(sn,i)
				plt.savefig(good_savedir + 'A_bayes_'+str(i)+'.png')
				plt.close()
		
		thres_size = plt_t.get_threshold_responses_size(sn)
		if thres_size>0:
			for i in range(thres_size):
				plt_t.plot_threshold_responses(sn,i)
				plt.savefig(good_savedir + 'B_threshold_' + str(i) + '.png')
				plt.close()
		
		all_size = plt_t.get_all_responses_size(sn)
		if all_size>0:
			for i in range(all_size):
				
				plt_t.plot_all_responses(sn,i)
				plt.savefig(sn_savedir+'B_all_'+str(i)+'.png')
				plt.close()

pool = Pool(2)     # Parallel computing
pool.map(run,list(t_sl))
pool.close()
pool.join()

