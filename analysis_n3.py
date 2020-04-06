import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from multiprocessing import Pool
import Data_analysis.file as myfile


databaselink = '/media/laojin/Elements/trigdata/'
yearlist = [2016]

def get_sample_dir_list(yearlist,databaselink):
	sample_dir_list = []
	for year in yearlist:
		topdir = databaselink + str(year) +'/'
		dirlist1 = os.listdir(topdir)
		for dirl in dirlist1:
			if os.path.isdir(topdir + dirl):
				sample_dir_list.append(topdir + dirl+'/')
	return sample_dir_list



def analysis_one(sample_link):
	'''
	
	:param sample_link:
	:return:
	'''
	NaI = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
	BGO = ['b0','b1']
	








