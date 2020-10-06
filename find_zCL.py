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
import Plume
imp.reload(Plume)
import graphics
imp.reload(graphics)

RunList = [i for i in config.tag if i not in config.exclude_bad]         #load a list of cases
runCnt = len(RunList)                                                  #count number of cases

#======================perform main analysis for all runs first===================
print('==========================LES analysis==========================')
#loop through all LES cases
all_plumes = []
for nCase,Case in enumerate(RunList):
    csdict = utils.prepCS(Case)      #load cross-wind integrated smoke and all other data

    #initalize a plume object
    plume = Plume.LESplume(Case)

    #get quasi-stationary profile
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= config.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells
    plume.get_zCL(pm)
    plume.classify()

    #estimate fire intensity
    plume.get_I(csdict['ghfx2D'],5000)


    #make plots, if necessary
    if config.plot_profiles:
        graphics.plot_profiles(plume,config.interpZ, csdict['w'][-1,:,:], csdict['temp'][-1,:,:])
    if config.plot_conservedvars:
        graphics.plot_conservedvars(plume,csdict['temp'][-1,:,:],pm)
    if config.plot_zcl:
        graphics.plot_zcl(plume,pm,csdict['ghfx'][-1,:],csdict['ghfx2D'][-1,:,:])

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

#obtain empirical parameter C
firstGuess = np.array(zPrimeGuess)[:,np.newaxis]
C, _, _, _ = np.linalg.lstsq(firstGuess, zPrimeTrue)
C = float(C)

#obtain bias correction factors
zCLmodel, zCLtrue = [],[]
for plume in penetrative_plumes:
    estimate = C*plume.Tau*plume.wf + plume.zs
    zCLmodel.append(estimate)
    zCLtrue.append(plume.zCL)

zCLmodel, zCLtrue = utils.plume_error(penetrative_plumes, C)
biasFit = linregress(zCLmodel,zCLtrue)

#plot model performance
graphics.injection_model(penetrative_plumes, C, biasFit)


#===========test iterative solution, do bias correction===============
imp.reload(Plume)
raw_error, unbiased_error, true_zCL = [], [], []
for plume in penetrative_plumes:
    raw_plume = Plume.MODplume(plume.name)
    raw_plume.I = plume.I
    raw_plume.iterate(C)
    unbiased_plume = Plume.MODplume(plume.name)
    unbiased_plume.I = plume.I
    unbiased_plume.iterate(C,biasFit)
    true_zCL.append(plume.zCL)
    raw_error.append(plume.zCL - raw_plume.zCL)
    unbiased_error.append(plume.zCL - unbiased_plume.zCL)

#plot bias correction statistics on iterative solution
graphics.bias_correction(raw_error, unbiased_error, true_zCL, figname='IterativeSolution')

#================explicit solution and dimensionless groups=========================
imp.reload(Plume)
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
    estimateH = C**(3/2.)*((plume.THs/(config.g*Gamma**3))**(1/4.)) * np.sqrt(plume.I/plume.zi**3)
    HStar.append(estimateH)                     #second dimensionless group
    cI.append(plume.I)              #save intensity for colorizing the plot

    #explicit solution
    raw_plume = Plume.MODplume(plume.name)
    raw_plume.I = plume.I
    raw_plume.explicit_solution(C, Gamma, ze)
    unbiased_plume = Plume.MODplume(plume.name)
    unbiased_plume.I = plume.I
    unbiased_plume.explicit_solution(C, Gamma, ze, biasFit)
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
    zCLmodel, zCLtrue = utils.plume_error(trainSubset, C)
    trialFit = linregress(zCLmodel, zCLtrue)
    Rstore.append(trialFit[2])

    test_error, test_truth, fuel_cat = [], [], []
    for nTest in TestSet:
        testPlume = Plume.MODplume(penetrative_plumes[nTest].name)
        testPlume.I = penetrative_plumes[nTest].I
        test_estimate = testPlume.iterate(C, trialFit)              #####STOPPPED HERE: NEED TO ADD POLYMORPHISM, TO PROVIDE OUTPUT TO FUNCITON (instead of assigning attributes)
        truth = penetrative_plumes[nTest].zCL
        test_truth.append(truth)
        test_error.append(truth - test_estimate)
        fuel_cat.append(utils.read_tag('F',np.array(penetrative_plumes[nTest].name)))

    error = test_truth  - test_error                         #calculate error between model and 'truth'
    ModelError.append(error)                                        #store model error
    TrueTrialZcl.append(test_truth)                           #store true subset
    TrialFuel.append(fuel_cat)

flatTrialZcl  = np.concatenate(TrueTrialZcl)                #flatten test array of injection heights
flatModelError = np.concatenate(ModelError)                     #flatten model error
flatTrialFuel = np.concatenate(TrialFuel)

#plot model sensitivity
plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(2, 2)
ax0 = plt.subplot(gs[0,0:])
plt.title('(a) TRIAL ERROR')
plt.boxplot(ModelError)
plt.hlines(0,0,11,colors='grey',linestyles='dashed')
ax0.set(xlabel='trial no.', ylabel='error in zCL [m]',ylim=[-200,200])

ax1 = plt.subplot(gs[1,0])
plt.title('(b) ALL TRIALS')
plt.boxplot(flatModelError)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel='all runs', ylabel='error in zCL [m]',ylim=[-200,200])

ax2 = plt.subplot(gs[1,1])
ax2.set(xlabel='R-value',ylabel='count' )
plt.title('(c) R-VALUE SENSITIVITY')
plt.hist(Rstore, bins=3)
ax2.set(xlabel='R-value',ylabel='count' )
plt.tight_layout()
plt.savefig(plume.figdir + 'injectionModel/ModelSensitivity.pdf')
plt.show()
# plt.close()


plt.figure()
plt.title('ERROR AS A FUNCTION OF FUEL (TRIALS)')
plt.scatter(flatTrialFuel,flatModelError)
plt.hlines(0,0,14,colors='grey',linestyles='dashed')
ax = plt.gca()
# for i, txt in enumerate(flatTrialName):
#     ax.annotate(txt, (flatTrialFuel[i], flatModelError[i]),fontsize=6)
ax.set(xlabel='fuel category', ylabel='error [m]',ylim=[-100,150])
plt.colorbar().set_label('$z_{CL} - z_s$ [m]')
plt.savefig(plume.figdir + 'injectionModel/FuelvsErrorHeight_TRIALS.pdf')
plt.show()
plt.close()

#
# #plot rearranged SOLUTION
# LHS = (zCL - zS) * (g*Omega/thetaS)
# RHS = (C**6) * ((Phi/(zi * Omega))**2)
# plt.figure()
# plt.title('REARRANGED FORM RS')
# plt.scatter(RHS, LHS)
# plt.gca().set(xlabel=r'$C^6\left(\frac{I}{z_i(\theta_{CL}-\theta_s)}\right)^2$',ylabel=r'$z\prime g\prime$')
# plt.savefig(plume.figdir + 'injectionModel/RearrangedGroupsRS.pdf')
# plt.show()
#
# #plot velocity comparison
# plt.figure()
# plt.title('COMPARE VELOCITIES')
# plt.scatter(wStarC/(Tau*C),(Phi/(zi * Omega)) )
# plt.gca().set(xlabel=r'$\widetilde{w_f}$',ylabel=r'$\frac{I}{z_i(\theta_{CL}-\theta_s)}$',aspect='equal',xlim = [0,20],ylim = [0,20])
# plt.savefig(plume.figdir + 'injectionModel/CompareWs.pdf')
# plt.show()
