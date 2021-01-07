
'''
nmoisseeva@eoas.ubc.ca
January 2018
This is sample code for interpolating and cross-sectioning synthetic plume data.
Note: this requires substantial memory and compute power.

'''

import numpy as np
from scipy.io import netcdf
from scipy import interpolate
import os.path
import wrf
import pickle
#====================INPUT===================
loc_tag = ['plume_tag1', 'plume_tag2', '...']


wrfdir = '../insert_path_to_raw_data/'
lvl = np.arange(0,3201,40)      #vertical levels to interpolate to
dx = 40.                        #horizontal grid spacing

cs = 10                         #+/- grids for cross-section
wi, wf= 25,375                  #window around the peak heat
ign_over = 20                    #number of history intervals to skip from beginning: 10min (600/15sec) or ignition (ceil(95sec / 15sec))
#=================end of input===============

print('INTERPOLATION AND AVERAGING SCRIPT FOR PLUME RUNS')
print('===================================')

for nCase,Case in enumerate(loc_tag):
    print('Examining case: %s ' %Case)
    #----------check for interpolated data----------------------------
    interppath = wrfdir + 'interp/wrfinterp_' + Case + '.npy'

    if os.path.isfile(interppath):

        print('Interpolated data found at: %s' %interppath)
        interpfile = open(interppath,'rb')
        interpdict = pickle.load(interpfile)   # load here the above pickle

    else:
        print('WARNING: no interpolated data found - generating: SLOW ROUTINE!')
        wrfpath = wrfdir + 'wrfout_' + Case
        wrfdata = netcdf.netcdf_file(wrfpath, mode ='r')
        ncdict = wrf.extract_vars(wrfdata, None, ('GRNHFX','W','QVAPOR','T','PHB','PH','U','P','PB','V','tr17_1'),meta=False)
        ncdict['PM25'] = ncdict.pop('tr17_1')

        #get height and destagger vars
        zstag = (ncdict['PHB'] + ncdict['PH'])//9.81
        z = wrf.destagger(zstag,1)
        u = wrf.destagger(ncdict['U'],3)
        w = wrf.destagger(ncdict['W'],1)
        v = wrf.destagger(ncdict['V'],2)

        #list variables to interpolate
        nT,nZ,nY,nX = np.shape(z)

        winterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
        uinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
        tinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
        vinterp = np.empty((nT,len(lvl),nY,nX)) * np.nan
        pminterp = np.empty((nT,len(lvl),nY,nX)) * np.nan

        for t in range(nT):
            print('.... tsetp = %s/%s' %(t,nT))
            for y in range(nY):
                for x in range(nX):
                    z_t = z[t,:,y,x]
                    ft = interpolate.interp1d(z_t,ncdict['T'][t,:,y,x],fill_value="extrapolate")
                    fw = interpolate.interp1d(z_t,w[t,:,y,x],fill_value="extrapolate")
                    fu = interpolate.interp1d(z_t,u[t,:,y,x],fill_value="extrapolate")
                    fv = interpolate.interp1d(z_t,v[t,:,y,x],fill_value="extrapolate")
                    fpm = interpolate.interp1d(z_t,ncdict['PM25'][t,:,y,x],fill_value="extrapolate")
                    winterp[t,:,y,x] = fw(lvl)
                    tinterp[t,:,y,x] = ft(lvl)
                    uinterp[t,:,y,x] = fu(lvl)
                    vinterp[t,:,y,x] = fv(lvl)
                    pminterp[t,:,y,x] = fpm(lvl)
                    interpdict = {'W':winterp, 'T':tinterp, 'U':uinterp,'V':vinterp,'PM25':pminterp,'GRNHFX': ncdict['GRNHFX']}
        writefile = open(interppath, 'wb')
        pickle.dump(interpdict, writefile, protocol=4)
        writefile.close()
        print('Interpolated data saved as: %s' %interppath)
        wrfdata.close()

    #convert and average data-----------------------------------
    ghfx = interpdict['GRNHFX']/1000.             #convert to kW
    temp = interpdict['T']+300.             #add perturbation and base temperature
    w = interpdict['W']
    u = interpdict['U']
    pm25 = interpdict['PM25']

    #get dimensions
    dimt, dimy, dimx = np.shape(ghfx)
    print('Dimensions of data: %s x %s x %s:'%(dimt,dimy,dimx))
    xsx = int(round(dimy/2.))

    var_list = ['ghfx','temp','w','u','pm25']
    csdict = {}

    for variable in var_list:
    #create fire cross-section
        if variable == 'ghfx':
            slab = vars()[variable][:,xsx-cs:xsx+cs,:]
            csdict[variable] = np.mean(slab,1)
        elif variable == 'pm25':
            slab = vars()[variable]
            csdict[variable] = np.sum(slab,2)   #tracers calculated as cross-wind intergrated totals
        else:
            slab = vars()[variable][:,:,xsx-cs:xsx+cs,:]
            csdict[variable] = np.mean(slab,2)
    csdict['ghfx2D'] = ghfx

    cspath = wrfdir + 'interp/wrfcs_' + Case + '.npy'
    np.save(cspath,csdict)
    print('Crosssection data saved as: %s' %cspath)
