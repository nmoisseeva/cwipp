# September 2020
#nmoisseeva@eoas.ubc.ca
#configuration settings

import numpy as np
import os.path
import glob

#====================INPUT BEGINS HERE================

#----------------------common settings----------------
#vertical level settings
zstep = 20                      #height interpolation step for analysis
interpZ = np.arange(0, 4001, zstep)     #set up interpolated vertical profile
BLfrac = 0.75                   #fraction of BL height to use as reference height z_s (default = 0.75)

#common analysis variables
g = 9.81                        #gravity constant
PMcutoff = 20                   #minimum PM value to define plume edge

#default model bias fit parameters
biasFit = [0.9195, 137.9193]

#figure storage directory
figdir = '/Users/nmoisseeva/code/pr_model/figs/'


#---------------------LES mode settings----------------

#path to interpolated cross-sectional plume data
wrfdir = '/Users/nmoisseeva/data/plume/main/interp/'

#grid spacing of interpolated LES data
dz = 40
dx = 40.
dy = 40.

#which runs to process
dirpath = wrfdir+'wrfcs_*'          #naming convention for cross-sections
exclude_runs = ['W5F4R0','W5F9R1','W5F8R3','W5F9R3','W5F1R3','W5F13R0','W5F1R7T','W5F8R7T','W5F9R7T']
exclude_bad = ['W5F4R0','W5F13R0']

#fireline configuration
ign_over = 20                   #number of history intervals exluded from analysis start (corresponding to ignition)

#which plots to make
plot_profiles = 1
plot_conservedvars = 1
plot_zcl = 1

#model cross-evaluation variables
trials = 10                     #number of times to rerun the model
testPortion = 0.2               #portion of data to reserve for testing the model

#------------------- do NOT modify unless sure----------
dirlist = glob.glob(dirpath)        #get directory listing
tag = [i[len(dirpath)-1:-4] for i in dirlist]
