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

# RunList=['W5F4R7T'] # for domain potting
RunList=['W5F12R6TE'] # for short vs tall plumes

# RunList = [i for i in config.tag if i not in config.exclude_bad]         #load a list of cases
runCnt = len(RunList)                                                  #count number of cases

#======================perform main analysis for all runs first===================
print('==========================LES analysis==========================')
#loop through all LES cases
all_plumes = []
for nCase,Case in enumerate(RunList):
    csdict = utils.prepCS(Case)      #load cross-wind integrated smoke and all other data

    #initalize a plume object
    plume = Plume.LESplume(Case)

    # #plot numerical setup
    # plt.figure(figsize=(7,3))
    # plt.imshow(csdict['ghfx2D'][15,:,:],cmap=plt.cm.YlOrRd, extent=[0,20000,0,10000],vmin=0, vmax = 300)
    # plt.gca().set(aspect='equal', xlabel='x (E/W) distance [m]', ylabel='y (N/S) distance [m]')
    # plt.colorbar(orientation='vertical', label=r'fire heat flux [$kW / m^2$]')
    # plt.tight_layout()
    # plt.savefig('domain_setup.pdf')
    # plt.show()


    pm = ma.masked_where(csdict['pm25'][-1,:,:] <= config.PMcutoff, csdict['pm25'][-1,:,:] ) #mask all non-plume cells

    #get quasi-stationary profile
    plume.get_zCL(pm, plot=config.plot_zcl, csdict=csdict)
    graphics.plot_conservedvars(plume,csdict['temp'][-1,:,:],pm)
