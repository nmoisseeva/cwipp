# March 2020
#nmoisseeva@eoas.ubc.ca
# This code partitions LES runs into model and test sets and applies the injection height parameterization
# Plotting shows model sensitivity and error distributions

import numpy as np
import matplotlib.pyplot as plt
# from scipy.io import netcdf
# import os.path
import imp
# from numpy import ma
# from matplotlib import ticker
# from scipy.signal import savgol_filter
# from scipy.stats import linregress
# from scipy.optimize import fsolve
# from matplotlib import gridspec
# from scipy.interpolate import interp1d
# import pickle

#import all common project variables
import config
imp.reload(config)
import utils
imp.reload(utils)
import cwipp
imp.reload(cwipp)
import graphics
imp.reload(graphics)



fires = np.load(config.input_data,allow_pickle=True).item()

for fire in fires:
	print('Processing fire: %s' %fire)
	case = cwipp.MODplume(fire)
	case.get_sounding(fires[case.name]['sounding'])
	case.I = fires[case.name]['I']
	case.iterate(config.biasFit)
	case.classify()
	if not case.penetrative:
		smoke_profile = '1'
	elif case.penetrative:

		'''
			determine plume top
			determine spread
			model full profile
		'''





#
# RunList = [i for i in config.tag if i not in config.exclude_bad]         #load a list of cases
# runCnt = len(RunList)                                                  #count number of cases
# # exclude_runs = ['W5F9R1','W5F8R3','W5F9R3','W5F1R3','W5F1R7T','W5F8R7T','W5F9R7T']
# #======================perform main analysis for all runs first===================
# print('==========================LES analysis==========================')
# #loop through all LES cases
# all_plumes = []
# for nCase,Case in enumerate(RunList):
#     csdict = utils.prepCS(Case)      #load cross-wind integrated smoke and all other data
#
#     #initalize a plume object
#     plume = cwipp.LESplume(Case)
#
#     #get quasi-stationary profile
#     pm = ma.masked_where(csdict['pm25'][-1,:,:] <= config.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells
#     plume.get_zCL(pm,plot=config.plot_zcl,csdict=csdict)
#     plume.classify()
#
#     #estimate fire intensity
#     plume.get_I(csdict['ghfx2D'],5000)
#
#
#     #make plots, if necessary
#     if config.plot_profiles:
#         graphics.plot_profiles(plume,config.interpZ, csdict['w'][-1,:,:], csdict['temp'][-1,:,:])
#     if config.plot_conservedvars:
#         graphics.plot_conservedvars(plume,csdict['temp'][-1,:,:],pm)
#
#     all_plumes.append(plume)
#
# #plot soundings
# graphics.plot_soundings(all_plumes)
#
# #pickle and save all plume data
# with open('plumes.pkl', 'wb') as f:
#     pickle.dump(all_plumes, f)
#
# # # Getting back the objects:
# # with open('plumes.pkl','rb') as f:
# #     all_plumes = pickle.load(f)
#
# #======================assess model performance========================
# print('Fitting all penetrative plumes using LES zCL')
# #make a list of penetrative plumes
# penetrative_plumes = []
# for plume in all_plumes:
#     if plume.penetrative:
#         plume.get_wf()
#         penetrative_plumes.append(plume)
#
# #get first estimate for the plume rise
# zPrimeGuess, zPrimeTrue  = [],[]
# for plume in penetrative_plumes:
#     zPrimeGuess.append(plume.Tau * plume.wf)
#     zPrimeTrue.append(plume.zCL - plume.zs)
#
# #obtain bias correction factors
# zCLmodel, zCLtrue = [],[]
# for plume in penetrative_plumes:
#     estimate = plume.Tau*plume.wf + plume.zs
#     zCLmodel.append(estimate)
#     zCLtrue.append(plume.zCL)
#
# zCLmodel, zCLtrue = utils.plume_error(penetrative_plumes)
# biasFit = linregress(zCLmodel,zCLtrue)
#
# #plot model performance
# graphics.injection_model(penetrative_plumes, biasFit)
#
# #===========test iterative solution, do bias correction===============
# raw_error, unbiased_error, true_zCL = [], [], []
# for plume in penetrative_plumes:
#     raw_plume = cwipp.MODplume(plume.name)
#     raw_plume.I = plume.I
#     raw_plume.iterate()
#     unbiased_plume = cwipp.MODplume(plume.name)
#     unbiased_plume.I = plume.I
#     unbiased_plume.iterate(biasFit)
#     true_zCL.append(plume.zCL)
#     raw_error.append(plume.zCL - raw_plume.zCL)
#     unbiased_error.append(plume.zCL - unbiased_plume.zCL)
# # for plume in all_plumes:
# #     if plume.name in exclude_runs:
# #         continue
# #     else:
# #         raw_plume = cwipp.MODplume(plume.name)
# #         raw_plume.I = plume.I
# #         raw_plume.iterate(C)
# #         unbiased_plume = cwipp.MODplume(plume.name)
# #         unbiased_plume.I = plume.I
# #         unbiased_plume.iterate(C,biasFit)
# #         true_zCL.append(plume.zCL)
# #         raw_error.append(plume.zCL - raw_plume.zCL)
# #         unbiased_error.append(plume.zCL - unbiased_plume.zCL)
#
# #plot bias correction statistics on iterative solution
# graphics.bias_correction(raw_error, unbiased_error, true_zCL, figname='IterativeSolution')
#
