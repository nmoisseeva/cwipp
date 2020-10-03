# September 2020
#nmoisseeva@eoas.ubc.ca
#graphings tools for plume rise analysis

import numpy as np
import matplotlib.pyplot as plt
# from scipy.io import netcdf
import os.path
# import imp
# from numpy import ma
from matplotlib import ticker
from scipy.stats import linregress
from scipy.interpolate import interp1d
import config
import utils
from matplotlib import gridspec


def plot_profiles(plume, interpZ, w, T):
    dimZ, dimX = np.shape(w)     #get shape of data
    metlvls = np.arange(0,dimZ*config.dz,config.dz)

    #get maximum vertical velocity
    Wmax = np.max(w,1)
    interpW = interp1d(metlvls, Wmax,fill_value='extrapolate')
    interpWmax = interpW(interpZ)

    #get maximum temperature
    Tmax = np.max(T,1)
    interpT = interp1d(metlvls, Tmax,fill_value='extrapolate')
    interpTmax = interpT(interpZ)

    #create a storage directory for plots
    figpath = config.figdir + 'profiles/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    #vertical concentration slice at donwind locations of wmax and qmax
    plt.figure(figsize=(10,4))
    plt.suptitle('%s' %plume.name)
    plt.subplot(121)
    ax1 = plt.gca()
    plt.title('MAX VERTICAL VELOCITY')
    plt.plot(interpWmax,interpZ,'k.-',label='$w_{max}$')
    plt.axhline(y = plume.zi, ls=':', c='darkgrey', label='zi')
    plt.axhline(y = plume.zCL,ls='--', c='red',label='z$_{CL}$')
    ax1.set(xlabel = 'velocity [m/s]', ylabel='height [m]',ylim = [0,metlvls[-1]] )
    plt.legend()

    plt.subplot(122)
    plt.title(r'PLUME vs AMBIENT $\theta$')
    ax2 = plt.gca()
    plt.plot(plume.sounding, interpZ, label=r'pre-ignition $\theta$ profile',c='lightblue')
    plt.plot(interpTmax,interpZ,c = 'orange',label=r'in-plume $\theta_{max}$')
    plt.axhline(y = plume.zi, ls=':', c='darkgrey', label='zi')
    plt.axhline(y = plume.zCL,ls='--', c='red',label='z$_{CL}$')
    ax2.set(xlabel = r'$\theta$ [K]', ylabel='height [m]' ,xlim = [285,330],ylim = [0,metlvls[-1]])
    plt.legend()
    plt.tight_layout()
    plt.savefig(figpath + 'profiles_%s.pdf' %plume.name)
    plt.close()

def plot_soundings(all_plumes):

    #create a storage directory for plots
    figpath = config.figdir + 'injection/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    plt.figure()
    plt.title('PRE-IGNITION ATMOSPHERIC PROFILES')
    leg_handles = []
    names_list = []
    soundings = []
    for plume in all_plumes:
        names_list.append(plume.name)
        soundings.append(plume.sounding)

    Rtag = np.array([i for i in utils.read_tag('R',names_list)])  #list of initialization rounds (different soundings)
    for R in set(Rtag):
        print(R)
        for Case in np.array(soundings)[Rtag==R]:
            lR = plt.plot(Case, config.interpZ, color='C%s' %R, linewidth=1, label='R%s' %R)
        leg_handles.extend(lR)
    plt.gca().set(xlabel='potential temperature [K]',ylabel='height [m]',xlim=[280,330],ylim=[0,2800])
    plt.legend(handles=leg_handles)
    plt.savefig(figpath + 'T0profiles.pdf')
    plt.close()


def injection_model(all_plumes,C,biasFit):

    #create a storage directory for plots
    figpath = config.figdir + 'injection/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    modelled_array = []
    true_array = []
    c_array = []
    for plume in all_plumes:
        modelled_zCL = C*plume.Tau*plume.wf + plume.zs
        modelled_array.append(modelled_zCL)
        true_array.append(plume.zCL)
        c_array.append(plume.I)

    plt.figure()
    plt.title('MODELLED SMOKE INJECTION HEIGHTS')
    ax = plt.gca()
    plt.scatter(modelled_array ,true_array,c=c_array,cmap =plt.cm.plasma)
    ax.set(ylabel = r'true $z_{CL}$ [m]', xlabel = r'model $z_{CL}$ [m]',aspect='equal')
    # ax.set(ylabel = r'true $z_{CL}$ [m]', xlabel = r'model $z_{CL}$ [m]',xlim = [400,3200], ylim = [400,3200])
    plt.colorbar(label=r'fireline intensity [K m$^2$/s]')
    plt.plot(np.sort(modelled_array),biasFit[0]*np.sort(modelled_array)+biasFit[1], color='black', label='linear regression fit')
    plt.plot(np.sort(modelled_array),np.sort(modelled_array), linestyle = 'dashed', color='grey', label='unity line')
    plt.legend()
    plt.savefig(figpath + 'MainInjection_allPlumes.pdf')
    plt.close()

def dimensionless_groups(HStar,zStar,cI):
    #create a storage directory for plots
    figpath = config.figdir + 'injection/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    plt.figure()
    plt.title('DIMENSIONLESS RELATIONSHIP')
    plt.scatter(HStar,zStar,c=cI,cmap=plt.cm.plasma)
    ax = plt.gca()
    # for i, txt in enumerate(RunList):
    #     ax.annotate(txt, (HStar[i], zStar[i]),fontsize=6)
    plt.gca().set(xlabel = r'$\overline{H}$', ylabel=r'$\overline{z}$')
    plt.colorbar(label=r'fireline intensity [Km$^2$/s]')
    plt.savefig(figpath + 'DimensionlessGroups.pdf')
    plt.close()

def bias_correction(raw_error, unbiased_error, true_zCL, figname='BiasCorrection'):

    #create a storage directory for plots
    figpath = config.figdir + 'injection/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    plt.figure(figsize=(9,6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
    ax0 = plt.subplot(gs[0])
    plt.title('(a) Error as f($z_{CL})$: RAW')
    plt.scatter(true_zCL,raw_error)
    plt.hlines(0,200,max(true_zCL)+100,colors='grey',linestyles='dashed')
    # for i, txt in enumerate(RunList):
    #     ax0.annotate(txt, (zCL[i], zCLerror[i]),fontsize=6)
    ax0.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
    ax1 = plt.subplot(gs[1])
    plt.title('(b) Error Statistics')
    plt.boxplot(raw_error)
    plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
    ax1.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
    ax2 = plt.subplot(gs[2])
    plt.title('(c) Error as f($z_{CL})$: BIAS CORRECTED')
    plt.scatter(true_zCL,unbiased_error)
    plt.hlines(0,200,max(true_zCL)+100,colors='grey',linestyles='dashed')
    ax2.set(xlabel=r'$z_{CL}$ [m]', ylabel='error [m]',ylim =[-350,350],xlim=[400,3200])
    ax3 = plt.subplot(gs[3])
    plt.title('(d) Error Statistics')
    plt.boxplot(unbiased_error)
    plt.hlines(0,0.5,1.5,colors='grey',linestyles='dashed')
    ax3.set(xlabel=r'$z_{CL}$',ylabel='error [m]',ylim = [-350,350], xticklabels=[''])
    plt.subplots_adjust(top=0.85)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(figpath + figname + '.pdf')
    plt.close()
