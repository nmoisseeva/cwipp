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
#loop through all LES cases
all_plumes = []
for nCase,Case in enumerate(RunList):
    csdict = utils.prepCS(Case)      #load cross-wind integrated smoke and all other data

    #initalize a plume object
    plume = Plume.LESplume(Case,config.interpZ)

    #get quasi-stationary profile
    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= config.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells
    plume.get_zCL(pm)
    plume.classify()

    #estimate fire intensity
    plume.get_I(csdict['ghfx2D'],5000)



    if config.plot_profiles:
        graphics.plot_profiles(plume,config.interpZ, csdict['w'][-1,:,:], csdict['temp'][-1,:,:])

    all_plumes.append(plume)


#pickle and save all plume data
with open('plumes.pkl', 'wb') as f:
    pickle.dump(all_plumes, f)


# Getting back the objects:
with open('plumes.pkl','rb') as f:
    all_plumes = pickle.load(f)

#======================assess model performance========================
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

#obtain emperical parameter C
firstGuess = np.array(zPrimeGuess)[:,np.newaxis]
C, _, _, _ = np.linalg.lstsq(firstGuess, zPrimeTrue)
C = float(C)

#obtain bias correction factors
zPrimeModel = C*modelGuess
zCLmodel, zCLtrue = [],[]
for plume in penetrative_plumes:
    estimate = C*plume.Tau*plume.wf + plume.zs
    zCLmodel.append(estimate)
    zCLtrue.append(plume.zCL)
biasFit = linregress(zCLmodel,zCLtrue)

# wStarC = Tau*((g*Phi*(zCL-zS))/(thetaS*zi))**(1/3.)
# firstGuessArray = wStarC[:,np.newaxis]
# C, _, _, _ = np.linalg.lstsq(firstGuessArray, zCL-zS)
# modelGuess = C*wStarC + zS
# biasFit = linregress(modelGuess,zCL)

#plot model performance
graphics.injection_model(penetrative_plumes, C, biasFit)

import graphics
imp.reload(graphics)

# plt.figure()
# plt.title('MODELLED SMOKE INJECTION HEIGHTS')
# ax = plt.gca()
# # plt.scatter(wStarC+zS,zCL,c=plume.read_tag('R',RunList),cmap =plt.cm.tab10)
# plt.scatter(modelGuess,zCL,c=Phi,cmap =plt.cm.plasma)
# # ax.set(ylabel = r'$z_{CL}$ [m]', xlabel = r'$C\tau_* \widetilde{w_f} + \frac{3}{4}z_i$ [m]',xlim = [400,3200], ylim = [400,3200])
# ax.set(ylabel = r'true $z_{CL}$ [m]', xlabel = r'injection model $z_{CL}$ [m]',xlim = [400,3200], ylim = [400,3200])
# plt.colorbar(label=r'fireline intensity [K m$^2$/s]')
# # plt.colorbar(label=r'atmospheric profile number')
# plt.plot(np.sort(modelGuess),biasFit[0]*np.sort(modelGuess)+biasFit[1], color='black', label='linear regression fit')
# plt.plot(np.sort(modelGuess),np.sort(modelGuess), linestyle = 'dashed', color='grey', label='unity line')
# plt.legend()
# plt.savefig(plume.figdir + 'injectionModel/NewInjectionTheory.pdf')
# plt.show()


plt.figure()
plt.title('PRE-IGNITION ATMOSPHERIC PROFILES')
leg_handles = []
Rtag = np.array([i for i in plume.read_tag('R',RunList)])  #list of initialization rounds (different soundings)
for R in set(Rtag):
    for Case in sounding[Rtag==R]:
        lR = plt.plot(Case, interpZ, color='C%s' %R, linewidth=1, label='R%s' %R)
    leg_handles.extend(lR)
plt.gca().set(xlabel='potential temperature [K]',ylabel='height [m]',xlim=[280,330],ylim=[0,2800])
plt.legend(handles=leg_handles)
plt.savefig(plume.figdir + 'T0profiles.pdf')
plt.show()

#================dimensionless fit=========================
R8idx = np.where(plume.read_tag('R', RunList)==8)[0][0]
exclude_idx = np.ndarray.flatten(np.arange(runCnt)!=R8idx)


Gamma = np.empty(len(RunList))* np.nan
thetaE =  np.empty(len(RunList))* np.nan
zE = np.empty(len(RunList))* np.nan
for nCase, Case in enumerate(RunList):
    ziidx = np.nanargmin(abs(interpZ - zi[nCase]))
    BLidx = int(ziidx*BLfrac)
    fittopidx = np.nanargmin(abs(interpZ - 2700))
    baselinedTheta = sounding[nCase,ziidx:fittopidx] - sounding[nCase,BLidx]
    GammaFit = linregress(interpZ[ziidx:fittopidx],baselinedTheta)
    Gamma[nCase] = GammaFit[0]
    thetaE[nCase] = sounding[nCase,BLidx]
    zE[nCase] = -GammaFit[1]/GammaFit[0]
zStar = (zCL - zE)/zi
# HStar = (3/4.)**(3/2.)*((thetaE/(g*Gamma**3))**(1/4.)) * np.sqrt((3/2.)*Phi/zi**3)
HStar = C**(3/2.)*((thetaE/(g*Gamma**3))**(1/4.)) * np.sqrt(Phi/zi**3)


plt.figure()
plt.title('DIMENSIONLESS RELATIONSHIP')
plt.scatter(HStar[exclude_idx],zStar[exclude_idx],c=Phi[exclude_idx],cmap=plt.cm.plasma)
ax = plt.gca()
# for i, txt in enumerate(RunList):
#     ax.annotate(txt, (HStar[i], zStar[i]),fontsize=6)
plt.gca().set(xlabel = r'$\overline{H}$', ylabel=r'$\overline{z}$')
plt.colorbar(label=r'fireline intensity [Km$^2$/s]')
plt.savefig(plume.figdir + 'injectionModel/DimensionlessGroups.pdf')
plt.show()


#===========iterative solution===============
zCLerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
zCLerrorBiased = np.empty((runCnt)) * np.nan          #parameterization error [m]
zCLcalc = np.empty((runCnt)) * np.nan          #parameterization error [m]

from scipy.optimize import root
for nCase,Case in enumerate(RunList):
    BLidx = np.nanargmin(abs(interpZ - BLfrac*zi[nCase]))


    # toSolveCase = lambda z : z - (C*BLfrac*zi[nCase] + wStarCFit[1]) - \
    #                 C * 1/(np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
    #                 (g*Phi[nCase]*(z-zi[nCase]*BLfrac)/(thetaS[nCase] * zi[nCase]))**(1/3.)

    toSolveCase = lambda z : z  - biasFit[1] - biasFit[0]*(zS[nCase] + \
                    C/(np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
                    (g*Phi[nCase]*(z-zi[nCase]*BLfrac)/(thetaS[nCase] * zi[nCase]))**(1/3.))

    #NOT bias-corrected
    # toSolveCaseBiased = lambda z : z - (BLfrac*zi[nCase])  - \
    #                 3./(4* np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
    #                 (g*Phi[nCase]*(z-zi[nCase]*BLfrac)*(3/2.)/(thetaS[nCase] * zi[nCase]))**(1/3.)
    #
    # toSolveCaseBiased = lambda z : z - (BLfrac*zi[nCase])  - \
    #                 C/(np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
    #                 (g*Phi[nCase]*(z-zi[nCase]*BLfrac)/(thetaS[nCase] * zi[nCase]))**(1/3.)

    toSolveCaseBiased = lambda z : z - zS[nCase]  - \
                    C/(np.sqrt(g*(sounding[nCase,int(z/zstep)] - thetaS[nCase])/(thetaS[nCase] * (z-zi[nCase]*BLfrac))))  * \
                    (g*Phi[nCase]*(z-zi[nCase]*BLfrac)/(thetaS[nCase] * zi[nCase]))**(1/3.)

    z_initial_guess = zi[nCase]                   #make initial guess BL height
    z_solution = fsolve(toSolveCase, z_initial_guess,factor=0.1)             #solve
    z_solutionBiased = fsolve(toSolveCaseBiased, z_initial_guess,factor=0.1)             #solve

    # zCLcalc[nCase] = (3/4)**6 * (thetaS[nCase]/g) * (((3/2.)* Phi[nCase]/zi[nCase])**2) * (thetaCL[nCase]-thetaS[nCase])**(-3) + BLfrac*zi[nCase]
    zCLerror[nCase] =  zCL[nCase]  - z_solution                              #store the solution
    zCLerrorBiased[nCase] =  zCL[nCase]  - z_solutionBiased                              #store the solution


plt.figure(figsize=(9,6))
plt.suptitle('ITERATIVE SOLUTION')
gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
ax0 = plt.subplot(gs[0])
plt.title('(a) Error as f($z_{CL})$: RAW')
plt.scatter(zCL,zCLerrorBiased)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax0.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax1 = plt.subplot(gs[1])
plt.title('(b) Error Statistics')
plt.boxplot(zCLerrorBiased)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
ax2 = plt.subplot(gs[2])
plt.title('(c) Error as f($z_{CL})$: BIAS CORRECTED')
plt.scatter(zCL,zCLerror)
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
ax2.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax3 = plt.subplot(gs[3])
plt.title('(d) Error Statistics')
plt.boxplot(zCLerror)
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax3.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
plt.show()
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'injectionModel/IterativeSolution.pdf')
plt.show()
plt.close()


# #========sensitivity test for zs values==========================
# plt.figure()
# ax = plt.gca()
# for BLfrac in np.arange(0.5, 0.78, 0.02):
#     omega = np.empty((len(RunList)))
#     thetas = np.empty((len(RunList)))
#     for nCase in range(len(RunList)):
#         sidx = np.nanargmin(abs(interpZ - BLfrac*zi[nCase]))
#         zCLidx = np.argmin(abs(interpZ - zCL[nCase]))
#         omega[nCase] = np.trapz(gradT0interp[nCase,sidx:zCLidx],dx=zstep)
#         thetas[nCase] = sounding[nCase,sidx]
#     tau = 1/ np.sqrt(g*omega/(thetas * (zCL-zS)))
#     wstar =  (3/4) *tau*(g*Phi*(zCL-zS)*(3/2.)/(thetas*zi))**(1/3.)       #attempt at zi_dependent threshold
#     wstarfit = linregress(wstar+zS,zCL)
#     print(wstarfit)
#     ax.scatter(BLfrac,wstarfit[2],c='C1')
# ax.set(xlabel='BL fraction', ylabel='R value of the fit')
# plt.tight_layout()
# plt.show()

#===========gamma solution===============
Gammaerror = np.empty((runCnt)) * np.nan          #parameterization error [m]
GammaerrorBiased = np.empty((runCnt)) * np.nan          #parameterization error [m]
for nCase,Case in enumerate(RunList):
    Gamma_solution = biasFit[0]*(C*((thetaE[nCase]/g)**(1/4.)) * ((Phi[nCase]/zi[nCase])**(0.5)) * (1/Gamma[nCase])**(3/4.) + zE[nCase]) + biasFit[1]
    Gammaerror[nCase] = zCL[nCase]  - Gamma_solution                               #store the solution
    Gamma_solutionBiased = (C**(3/2.)*(thetaE[nCase]/g)**(1/4.)) * ((Phi[nCase]/zi[nCase])**(0.5)) * (1/Gamma[nCase])**(3/4.) +zE[nCase]
    GammaerrorBiased[nCase] =  zCL[nCase]  - Gamma_solutionBiased                              #store the solution

plt.figure(figsize=(9,6))
plt.suptitle('EXPLICIT SOLUTION')
gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
ax0 = plt.subplot(gs[0])
plt.title(r'(a) Error as f($z_{CL})$: RAW')
plt.scatter(zCL[exclude_idx],GammaerrorBiased[exclude_idx])
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
# for i, txt in enumerate(RunList):
#     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
ax0.set(xlabel=r'$z_{CL}$ [m] ', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax1 = plt.subplot(gs[1])
plt.title('(b) Error Statistics')
plt.boxplot(GammaerrorBiased[exclude_idx])
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax1.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
ax2 = plt.subplot(gs[2])
plt.title(r'(c) Error as f($z_{CL})$: BIAS CORRECTED')
plt.scatter(zCL[exclude_idx],Gammaerror[exclude_idx])
plt.hlines(0,200,3200,colors='grey',linestyles='dashed')
ax2.set(xlabel=r'$z_{CL}$ [m] ', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
ax3 = plt.subplot(gs[3])
plt.title('(d) Error Statistics')
plt.boxplot(Gammaerror[exclude_idx])
plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
ax3.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
plt.subplots_adjust(top=0.85)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(plume.figdir + 'injectionModel/ExplicitSolution.pdf')
plt.show()
plt.close()



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


#plot rearranged SOLUTION
LHS = (zCL - zS) * (g*Omega/thetaS)
RHS = (C**6) * ((Phi/(zi * Omega))**2)
plt.figure()
plt.title('REARRANGED FORM RS')
plt.scatter(RHS, LHS)
plt.gca().set(xlabel=r'$C^6\left(\frac{I}{z_i(\theta_{CL}-\theta_s)}\right)^2$',ylabel=r'$z\prime g\prime$')
plt.savefig(plume.figdir + 'injectionModel/RearrangedGroupsRS.pdf')
plt.show()

#plot rearranged SOLUTION-with Tau
zCLrearranged = C**6 * (thetaS/g) * (Phi/zi)**2 * Omega**(-3.)+ zS
plt.figure()
plt.title('REARRANGED FORM NM')
plt.scatter(zCLrearranged, zCL)
plt.gca().set(xlabel=r'$C^6\left[\frac{\theta_s}{g}\right] \left[\frac{I}{z_i}\right]^2 \left[\theta_{CL}-\theta_s\right]+z_s$',ylabel=r'$z_{CL}$')
plt.savefig(plume.figdir + 'injectionModel/RearrangedGroupsNM.pdf')
plt.show()

#plot velocity comparison
plt.figure()
plt.title('COMPARE VELOCITIES')
plt.scatter(wStarC/(Tau*C),(Phi/(zi * Omega)) )
plt.gca().set(xlabel=r'$\widetilde{w_f}$',ylabel=r'$\frac{I}{z_i(\theta_{CL}-\theta_s)}$',aspect='equal',xlim = [0,20],ylim = [0,20])
plt.savefig(plume.figdir + 'injectionModel/CompareWs.pdf')
plt.show()
