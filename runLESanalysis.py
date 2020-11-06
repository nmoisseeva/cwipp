# March 2020
#nmoisseeva@eoas.ubc.ca
# This code partitions LES runs into model and test sets and applies the injection height parameterization
# Plotting shows model sensitivity and error distributions

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os.path
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.optimize import fsolve
from matplotlib import gridspec
from scipy.interpolate import interp1d
import pickle

#import all common project variables
import config
imp.reload(config)
import utils
imp.reload(utils)
import cwipp
imp.reload(cwipp)
import graphics
imp.reload(graphics)

RunList = [i for i in config.tag if i not in config.exclude_bad]         #load a list of cases
runCnt = len(RunList)                                                  #count number of cases
# exclude_runs = ['W5F9R1','W5F8R3','W5F9R3','W5F1R3','W5F1R7T','W5F8R7T','W5F9R7T']
#======================perform main analysis for all runs first===================
print('==========================LES analysis==========================')
#loop through all LES cases
all_plumes = []
for nCase,Case in enumerate(RunList):
    csdict = utils.prepCS(Case)      #load cross-wind integrated smoke and all other data

    #initalize a plume object
    plume = cwipp.LESplume(Case)

    #calculate sounding-related variables
    T0 = np.load(config.wrfdir + 'profiles/profT0' + name + '.npy')
    plume.get_sounding(T0) #!!!!! NEED TO TEST THIS

    #get quasi-stationary profile
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= config.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells
    plume.get_zCL(pm,plot=config.plot_zcl,csdict=csdict)
    plume.classify()

    #estimate fire intensity
    plume.get_I(csdict['ghfx2D'],5000)


    #make plots, if necessary
    if config.plot_profiles:
        graphics.plot_profiles(plume,config.interpZ, csdict['w'][-1,:,:], csdict['temp'][-1,:,:])
    if config.plot_conservedvars:
        graphics.plot_conservedvars(plume,csdict['temp'][-1,:,:],pm)

    all_plumes.append(plume)

#plot soundings
graphics.plot_soundings(all_plumes)

#pickle and save all plume data
with open('plumes.pkl', 'wb') as f:
    pickle.dump(all_plumes, f)

# # Getting back the objects:
# with open('plumes.pkl','rb') as f:
#     all_plumes = pickle.load(f)

#======================assess model performance========================
print('Fitting all penetrative plumes using LES zCL')
#make a list of penetrative plumes
penetrative_plumes = []
for plume in all_plumes:
    if plume.penetrative:
        plume.get_wf()
        penetrative_plumes.append(plume)

#get first estimate for the plume rise
zPrimeGuess, zPrimeTrue  = [],[]
for plume in penetrative_plumes:
    zPrimeGuess.append(plume.Tau * plume.wf)
    zPrimeTrue.append(plume.zCL - plume.zs)

# #obtain empirical parameter C
# firstGuess = np.array(zPrimeGuess)[:,np.newaxis]
# C, _, _, _ = np.linalg.lstsq(firstGuess, zPrimeTrue)
# C = float(C)
# C = 1

#obtain bias correction factors
zCLmodel, zCLtrue = [],[]
for plume in penetrative_plumes:
    estimate = plume.Tau*plume.wf + plume.zs
    zCLmodel.append(estimate)
    zCLtrue.append(plume.zCL)

zCLmodel, zCLtrue = utils.plume_error(penetrative_plumes)
biasFit = linregress(zCLmodel,zCLtrue)

#plot model performance
graphics.injection_model(penetrative_plumes, biasFit)

#===========test iterative solution, do bias correction===============
raw_error, unbiased_error, true_zCL = [], [], []
for plume in penetrative_plumes:
    raw_plume = cwipp.MODplume(plume.name)
    raw_plume.I = plume.I
    raw_plume.iterate()
    unbiased_plume = cwipp.MODplume(plume.name)
    unbiased_plume.I = plume.I
    unbiased_plume.iterate(biasFit)
    true_zCL.append(plume.zCL)
    raw_error.append(plume.zCL - raw_plume.zCL)
    unbiased_error.append(plume.zCL - unbiased_plume.zCL)
# for plume in all_plumes:
#     if plume.name in exclude_runs:
#         continue
#     else:
#         raw_plume = cwipp.MODplume(plume.name)
#         raw_plume.I = plume.I
#         raw_plume.iterate(C)
#         unbiased_plume = cwipp.MODplume(plume.name)
#         unbiased_plume.I = plume.I
#         unbiased_plume.iterate(C,biasFit)
#         true_zCL.append(plume.zCL)
#         raw_error.append(plume.zCL - raw_plume.zCL)
#         unbiased_error.append(plume.zCL - unbiased_plume.zCL)

#plot bias correction statistics on iterative solution
graphics.bias_correction(raw_error, unbiased_error, true_zCL, figname='IterativeSolution')

#================explicit solution and dimensionless groups=========================
print('Fitting dimensionless groups')
zStar, HStar, cI = [], [], []
raw_error, unbiased_error, true_zCL = [], [], []
for plume in penetrative_plumes:
    if utils.read_tag('R',[plume.name])==8:
        print('..... skipping %s' %plume.name)
        continue
    #dimenionsless groups
    i_zi = np.nanargmin(abs(config.interpZ - plume.zi))
    i_zs = np.nanargmin(abs(config.interpZ - plume.zs))
    i_top = np.nanargmin(abs(config.interpZ - (plume.zi+1000)))     #estimate lapse rate over a km above zi
    baselinedTH = plume.sounding[i_zi:i_top] - plume.sounding[i_zs]
    GammaFit = linregress(config.interpZ[i_zi:i_top],baselinedTH)
    Gamma = GammaFit[0]
    ze = -GammaFit[1]/GammaFit[0]
    zStar.append((plume.zCL - ze)/plume.zi)     #first dimensionless group
    estimateH = ((plume.THs/(config.g*Gamma**3))**(1/4.)) * np.sqrt(plume.I/plume.zi**3)
    HStar.append(estimateH)                     #second dimensionless group
    cI.append(plume.I)              #save intensity for colorizing the plot

    #explicit solution
    raw_plume = cwipp.MODplume(plume.name)
    raw_plume.I = plume.I
    raw_plume.explicit_solution(Gamma, ze)
    unbiased_plume = cwipp.MODplume(plume.name)
    unbiased_plume.I = plume.I
    unbiased_plume.explicit_solution(Gamma, ze, biasFit)
    true_zCL.append(plume.zCL)
    raw_error.append(plume.zCL - raw_plume.zCL)
    unbiased_error.append(plume.zCL - unbiased_plume.zCL)


#plot dimensionless performance
graphics.dimensionless_groups(HStar,zStar,cI)

#plot bias correction statistics on explicit solution
graphics.bias_correction(raw_error, unbiased_error, true_zCL, figname='ExplicitSolution')

#======================train and test regression model===================

#create storage arrays for R values, modelled zCL, model error and trail subsets of true zCL derived from data
ModelError, TrialZcl, TrialFuel, Rstore = [], [], [], []

for nTrial in range(config.trials):
    #split runs into train and test datasets
    RandomSample = np.random.binomial(1,config.testPortion,len(penetrative_plumes))                 #pick a random set of approximately 20% of data
    TrainSet = np.where(RandomSample==1)[0]
    TestSet = np.where(RandomSample==0)[0]

    #linear regression using training data subset only
    trainSubset = []
    for i in TrainSet:
        trainSubset.append(penetrative_plumes[i])
    zCLmodel, zCLtrue = utils.plume_error(trainSubset)
    trialFit = linregress(zCLmodel, zCLtrue)
    Rstore.append(trialFit[2])

    test_error, test_truth, fuel_cat = [], [], []
    for nTest in TestSet:
        testPlume = cwipp.MODplume(penetrative_plumes[nTest].name)
        testPlume.I = penetrative_plumes[nTest].I
        test_estimate, TH_dump = testPlume.iterate(trialFit, argout=True)              #####STOPPPED HERE: NEED TO ADD POLYMORPHISM, TO PROVIDE OUTPUT TO FUNCITON (instead of assigning attributes)
        truth = penetrative_plumes[nTest].zCL
        test_truth.append(truth)
        test_error.append(truth - test_estimate)
        fuel_cat.append(utils.read_tag('F',[penetrative_plumes[nTest].name])[0])

    # error = np.array(test_truth)  - np.array(test_error)                         #calculate error between model and 'truth'
    ModelError.append(test_error)                              #store model error
    TrialZcl.append(test_truth)                           #store true subset
    TrialFuel.append(fuel_cat)

flatTrialZcl  = np.concatenate(TrialZcl)                #flatten test array of injection heights
flatModelError = np.concatenate(ModelError)                     #flatten model error
flatTrialFuel = np.concatenate(TrialFuel)

graphics.model_sensitivity(ModelError,flatModelError,Rstore)
graphics.fuel_error(flatTrialFuel,flatModelError,flatTrialZcl)
