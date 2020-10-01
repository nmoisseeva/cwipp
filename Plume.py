import numpy as np
from config import *
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

class Plume:
    """
    Plume class

    ...

    Attributes
    ----------
    name : str
        plume name
    zi : float
        boundary layer height
    zs : float
        refernce height (zi * BLfrac)
    profile : ndarray
        1D vector corresponding to quasi-stationary downwind PM profile
    quartiles: ndarray
        2D array with columns corresponding to Q1 and Q3 profiles


    Methods
    -------
    get_profile(pm):
        Finds quasi-stationary downwind profile and its IQR
    """

    def __init__(self, name, zi, zs):
        """
        Constructs the plume object with some inital attributes

        Parameters
        ----------
            name : str
                plume name
            zi : float
                boundary layer height
            zs : float
                refernce height (zi * BLfrac)
        """
        self.name = name
        self.zi = zi
        self.zs = zs

    def get_profile(self, pm):
        """
        Find quasi-stationary profile

        Parameters
        ----------
            pm : ndarray
                2D array (z,x) of pm cross-section
        """

        #set up dimensions
        dimZ, dimX = np.shape(pm)     #get shape of data
        pmlvls = np.arange(0,dimZ,dz)

        #locate centerline
        ctrZidx = pm.argmax(0)                          #locate maxima along height
        ctrXidx = pm.argmax(1)                          #locate maxima downwind
        pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])                #get concentration along the centerline

        xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)                       #get location of maximum centerline height
        centerline = ma.masked_where(pmlvls[ctrZidx] == 0, pmlvls[ctrZidx])         #make sure centerline is only calculated inside the plume
        centerline.mask[:int(1000/dx)] = True

        filter_window = max(int(utlis.read_tag('W',[Plume.name])*10+1),51)
        smoothCenterline = savgol_filter(centerline, filter_window, 3)              #smooth centerline height

        #calculate concentration changes along the centerline
        dPMdX = pmCtr[1:]-pmCtr[0:-1]
        smoothPM = savgol_filter(dPMdX, filter_window, 3)

        #find where profile is quasi-stationary
        stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                                nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                nX > np.nanargmax(smoothPM) and\
                                nX > np.nanargmax(centerline) +10 and\
                                centerline[nX] < pmlvls[-1]-200 and \
                                nX > np.nanargmax(smoothCenterline)+10 else \
                                False for nX in range(dimX-1) ]
        if sum(stablePMmask) == 0:
            stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                    abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                                    nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                    nX > np.nanargmax(centerline) +10 and\
                                    nX > np.nanargmax(smoothPM) else\
                                    False for nX in range(dimX-1) ]
        if sum(stablePMmask) == 0:
            stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                    abs(smoothCenterline[nX+1]-smoothCenterline[nX]) < 5 and \
                                    nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                    nX > np.nanargmax(smoothPM) else\
                                    False for nX in range(dimX-1) ]

        stablePM = pm[:,1:][:,stablePMmask]
        stableProfile = np.median(stablePM,1)

        #find IQR
        pmQ1 = np.percentile(stablePM,25,axis = 1)
        pmQ3 = np.percentile(stablePM,75,axis = 1)
        interpQ1 = interp1d(pmlvls,pmQ1,fill_value='extrapolate')(interpZ)
        interpQ3 = interp1d(pmlvls,pmQ3,fill_value='extrapolate')(interpZ)

        #save attributes
        self.profile = interp1d(pmlvls,stableProfile,fill_value='extrapolate')(interpZ)
        self.quartiles = np.array([Q1,Q3])
