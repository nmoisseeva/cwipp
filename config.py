# September 2020
#nmoisseeva@eoas.ubc.ca
#configuration settings

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import imp
from numpy import ma
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import os.path
import glob


#paths
wrfdir = '/Users/nmoisseeva/data/plume/main/interp/'
figdir = '/Users/nmoisseeva/code/plume/main/clean_code/figs/'

#grid spacing of interpolated LES data
dz = 40
dx = 40.
dy = 40.

#which runs to process
dirpath = wrfdir+'interp/wrfcs_*'
dirlist = glob.glob(dirpath)
tag = [i[len(dirpath)-1:-4] for i in dirlist]
exclude_runs = ['W5F4R0','W5F9R1','W5F8R3','W5F9R3','W5F1R3','W5F13R0','W5F1R7T','W5F8R7T','W5F9R7T']
fireline_runs = ['W4F7R4L1','W4F7R4','W4F7R4L4']
exclude_bad = ['W5F4R0','W5F13R0']


#common analysis variables
g = 9.81                        #gravity constant
PMcutoff = 20                   #minimum value for plume edge
zstep = 20                      #height interpolation step for analysis
interpZ = np.arange(0, 4001, zstep) #set up interpolated vertical profile
BLfrac = 0.75                   #fraction of BL height to use as reference height z_s (default = 0.75)




# lvl = np.arange(0,2800,dz)	 	#vertical levels in m
# lvltall = np.arange(0,3201,dz)  #vertical levels for *T runs

#fireline configuration
ign_over = 20                   #number of history intervals exluded from analysis start (corresponding to ignition)

# cs = 10                         #+/- grids for cross-section
# wi, wf = 25, 375
# fireline = 4000.,6000.          #fireline start and end in meters



#model evaluation
rxcadredata = '/Users/nmoisseeva/data/plume/rxcadre/wrfout_test_main'
rxdispersion = '/Users/nmoisseeva/data/rxcadre/dispersion/RDS-2014-0015/Data/SmokeDispersion_L2G_20121110.csv'
rxsfchgt = 62               #surface height MSL (m)
rxsounding = '/Users/nmoisseeva/code/plume/main/input_sounding_rxcadre'
rxinterpCO2 = '/Users/nmoisseeva/data/plume/main/CO2_interp_1250-1300pm.npy'
rxlvl = np.arange(0,1700,20)
