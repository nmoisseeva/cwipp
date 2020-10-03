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
Rstore = np.empty((trials)) * np.nan
ModelError = []
TrueTrialZcl = []
TrialFuel = []
TrialName = []

for nTrial in range(trials):
    #split runs into train and test datasets
    TestFlag = np.random.binomial(1,testPortion,runCnt)                 #pick a random set of approximately 20% of data
    testCnt = sum(TestFlag)         #count how many runs ended up as test dataset

    #linear regression using training data subset only
    trialFit = linregress(C*wStarC[TestFlag==0]+zS[TestFlag==0],zCL[TestFlag==0])
    print('Sum of residuals using TRAINING data: %0.2f' %trialFit[2])
    Rstore[nTrial] = trialFit[2]        #store trial value

    #plot individual trial results
    fig = plt.figure()
    plt.suptitle('REGRESSION MODEL: ALL [R=%0.2f] vs TRAIN DATA [R=%0.2f]' %(biasFit[2], trialFit[2]))
    ax=plt.gca()
    plt.scatter(C*wStarC[TestFlag==0]+zS[TestFlag==0], zCL[TestFlag==0], c='C2', label='training data')
    plt.scatter(C*wStarC[TestFlag==1]+zS[TestFlag==1], zCL[TestFlag==1], c='C1', label='test data')
    plt.plot(np.sort(C*wStarC+zS),biasFit[0]* np.sort(C*wStarC+zS) + biasFit[1],  c='grey', label='all data')
    plt.plot(np.sort(C*wStarC+zS),trialFit[0]* np.sort(C*wStarC+zS) + trialFit[1], c='C2', label='training data regression fit')

    ax.set(xlabel=r'$\frac{3}{4}\tau \widetilde{w_{f*}} + \frac{3}{4}z_i$ [m/s]',ylabel='zCL [m]')
    plt.legend()
    plt.savefig(plume.figdir + 'injectionModel/trials/ALLvsTRAINdataTrial%s.pdf' %nTrial )
    plt.show()
    plt.close()

    #solve numerically for each run belonging to the Test subgroup
    zCLmodel = np.empty((testCnt)) * np.nan

    for nTest in range(testCnt):
        testIdx = np.where(TestFlag==1)[0][nTest]                   #get index of test run in the LES subset

        # toSolve = lambda z : z - (trialFit[0]*BLfrac*zi[testIdx] + trialFit[1]) - \
        #             trialFit[0] * 1/(np.sqrt(g*(sounding[testIdx,int(z/zstep)] - thetaS[testIdx])/(thetaS[testIdx] * (z-zi[testIdx]*BLfrac))))  * \
        #             (g*Phi[testIdx]*(z-zi[testIdx]*BLfrac)/(thetaS[testIdx] * zi[testIdx]))**(1/3.)

        toSolve = lambda z : z - trialFit[1] - trialFit[0]*(zS[testIdx] + \
                    C/(np.sqrt(g*(sounding[testIdx,int(z/zstep)] - thetaS[testIdx])/(thetaS[testIdx] * (z-zi[testIdx]*BLfrac))))  * \
                    (g*Phi[testIdx]*(z-zi[testIdx]*BLfrac)/(thetaS[testIdx] * zi[testIdx]))**(1/3.))

        z_initial_guess = zi[testIdx]                    #make initial guess BL height
        z_solution = fsolve(toSolve, z_initial_guess,factor=0.1)               #solve

        zCLmodel[nTest] = z_solution                                #store the solution
        print('%s solution is zCL = %0.2f' % (np.array(RunList)[testIdx],z_solution))
        print('...True value: %0.2f ' %zCL[testIdx])
    error = zCL[TestFlag==1]  - zCLmodel                          #calculate error between model and 'truth'
    ModelError.append(error)                                        #store model error
    TrueTrialZcl.append(zCL[TestFlag==1])                           #store true subset
    category = plume.read_tag('F',np.array(RunList)[TestFlag==1])
    TrialFuel.append(category)
    TrialName.append(np.array(RunList)[TestFlag==1])
#======================plot model sensitivity===================

flatTrueTrialZcl  = np.concatenate(TrueTrialZcl)                #flatten test array of injection heights
flatModelError = np.concatenate(ModelError)                     #flatten model error
flatTrialFuel = np.concatenate(TrialFuel)
flatTrialName = np.concatenate(TrialName)


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
