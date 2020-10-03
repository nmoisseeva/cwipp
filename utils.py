# September 2020
#nmoisseeva@eoas.ubc.ca
#fucntions and utilities

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
from scipy.interpolate import interp1d
import re
import config

# find zi
def get_zi(T0,dz):
    r'''
    Retrieve the height of boundary layer top $z_i$, based on gradient of potential temperature lapse rate.
    ...

    Parameters:
    -----------
    T0 : ndarray
        1D array representing potential temperature sounding [K]
    dz : float
        vertical grid spacing [m]

    Returns:
    --------
    $z_i$ : float
        boundary layer height [m]
    '''
    dT = T0[1:]-T0[0:-1]
    gradT = dT[1:] - dT[0:-1]
    surface = round(200/dz)
    zi_idx = np.argmax(gradT[surface:]) + surface                 #vertical level index of BL top
    zi = dz * zi_idx
    return zi


#tag reading function read_tag(variable type, string array)
def read_tag(condition, RunList):
    r'''
    Extract condition tags from all runs
    ...
    Parameters:
    -----------
    condition : str
        letter corresponding to test condition ('W'/'F'/'R'/'L')
    RunList : list
        list of strings containing plume tags

    Returns:
    --------
    out_array : ndarray
        1D array containing condition values for RunList
    '''
    out_array = []
    for nTag,tag in enumerate(RunList):
        letters = re.split('\d+', tag)
        numbers = re.findall('\d+', tag)
        if condition=='W':
            out_array.append(int(numbers[0]))
        elif condition=='F':
            out_array.append(int(numbers[1]))
        elif condition=='R':
            out_array.append(int(numbers[2]))
        elif condition=='L':
            if letters[-2]=='L':        #-2 accounts for empty string at the end of the file
                out_array.append(int(numbers[3]))
            else:
                out_array.append(2)
    out_array = np.array(out_array)
    return out_array


#load crosssection dictionary, extract profiles them and load them in to variables
def prepCS(Case):
    r'''
    Load simulation cross-section; prepare and save pre-ignition potential temperature and horizontal velocity profiles.

    Parameters:
    -----------
    Case : str
        plume name; assumed naming convention for cross-sectional data: wrfcs_Case.npy, \\
        containing dictionary keys 'temp' and 'u' for potentail temperature and horizontal velocity, respectively.

    Returns:
    --------
    csdict : dict
        dictionary containing cross-sectional plume data
    '''
    #load cross section
    cspath = config.wrfdir + 'wrfcs_' + Case + '.npy'
    print('Opening data: %s' %cspath)
    csdict = np.load(cspath, allow_pickle=True).item()

    #create a subfolder in the data folder to store profiles
    profilespath = config.wrfdir + 'profiles/'
    if not os.path.exists(profilespath):
        os.makedirs(profilespath)

    #save initial profiles
    profpathT = profilespath + 'profT0' + Case + '.npy'
    profpathU = profilespath + 'profU0' + Case + '.npy'
    if not os.path.isfile(profpathT) or not os.path.isfile(profpathU):
        profileT = np.mean(csdict['temp'][0,:,:],1)
        np.save(profpathT,profileT)

        profileU = np.mean(csdict['u'][0,:,:],1)
        np.save(profpathU,profileU)
        print('...Generated new profiles' )
    return csdict
