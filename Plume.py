import numpy as np
import config
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import utils


class Plume:
    """
    Parent Plume class

    ...

    Attributes
    ----------
    name : str
        plume name
    zi : float
        boundary layer height [m]
    zs : float
        refernce height (zi * BLfrac) [m]
    sounding: ndarray
        vertical potential temperature sounding on interpolated analysis levels [K]
    THs : float
        ambient potential tempreature at reference height zs [K]
    I : float
        fireline intensity parameter [K m2 s-1]

    Methods:
    -------
    get_I(flux2D, *Us):
        Finds cross-wind fireline intensity parameter I
    """

    def __init__(self, name, interpZ):
        """
        Constructs the plume object with some inital attributes
        ...

        Parameters:
        -----------
        name: str
            plume name
        interpZ: ndarray
            1D array containing levels AGL [m] on which to perform analysis
        """

        #get initial raw sounding (from cross-section wrfcs data)
        T0 = np.load(config.wrfdir + 'profiles/profT0' + name + '.npy')    #load initial temperature profile

        #get BL height
        zi = utils.get_zi(T0,config.dz)                    #calculate BL height
        zs = zi * BLfrac

        #interpolate sounding to analysis levels
        metlvls = np.arange(0,len(T0)*config.dz,config.dz)
        interpT= interp1d(metlvls,T0,fill_value='extrapolate')
        T0interp = interpT(interpZ)
        i_zs = np.argmin(abs(interpZ - zs))
        THs = T0interp[i_zs]

        self.name = name
        self.zi = zi
        self.zs = zs
        self.sounding = T0interp
        self.THs = THs


    def get_I(self, flux2D, depth, *Us):
        """
        Finds cross-wind fireline intensity parameter I
        ...
        Parameters:
        -----------
        flux2D: ndarray
            3D (time,y,x) or 2D (y,x) array containing heat flux values [kW m-2]
        depth: float
            maximum cross-wind depth of the fire over the entire timespan [m]
        Us: float, optional
            surface wind direction [deg, relative to y axis] NOT CURRENTLY IMPLEMENTED!

        Returns:
        -------
        I : float
            fireline intensity parameter [K m2 s-1]
        """

        #confirm there are sufficiant dimensions
        dims = np.shape(flux2D)
        if len(dims) > 3:
            raise ValueError('Too many dimensions for heat flux data')
        elif len(dims)<3:
            raise ValueError('Too few dimensions: need 3D array (time,y,x)')

        #mask and pad the heat source ------------------------
        upwind_padding = int(depth/config.dx)
        downwind_padding = int(1000/config.dx)              #assumes ground is not heated beyont 1km downwind
        masked_flux = ma.masked_less_equal(np.pad(flux2D,((0,0),(0,0),(upwind_padding,0)), 'constant',constant_values=0),1)

        cs_flux = np.nanmean(masked_flux,1)                         #get mean cross section for each timestep
        fire = []                                                   #create storage arrage
        fxmax = np.argmax(cs_flux,axis=1)                           #get location of max heat for each timestep
        for nP, pt in enumerate(fxmax[config.ign_over:]):            #excludes steps containing ignition
            subset = cs_flux[config.ign_over+nP,pt-upwind_padding:pt+downwind_padding]     #set averaging window around a maximum
            fire.append(subset)

        meanFire = np.nanmean(fire,0)                               #calculate mean fire cross section
        ignited = np.array([i for i in meanFire if i > 0.5])        #consider only cells that have heat flux about 500 W/m2
        I = np.trapz(ignited, dx = config.dx) * 1000 / ( 1.2 * 1005)    #calculate Phi by integrating kinematic heat flux along x (Km2/s)

        self.I = I

class LESplume(Plume):
    """
    Child Plume class used for operating on simulated plumes with full fields available (i.e. non-predictive mode)
    ...
    Attributes
    ----------
    profile : ndarray
        1D vector corresponding to quasi-stationary downwind PM profile [concentration]
    quartiles : ndarray
        2D array with columns corresponding to Q1 and Q3 profiles [concentration]
    zCL : float
        plume injection height [m]
    THzCL : float
        ambient potential temperature at zCL [K]

    Methods
    -------
    get_zCL(pm):
        Finds quasi-stationary downwind profile and its IQR, extracts injection height and associated variables
    """

    def get_zCL(self, pm):
        """
        Finds quasi-stationary downwind profile and its IQR, extracts injection height and associated variables

        Parameters
        ----------
        pm : ndarray
            2D array (z,x) of pm cross-section

        Returns:
        --------
        profile : ndarray
            1D vector corresponding to quasi-stationary downwind PM profile
        quartiles: ndarray
            2D array with columns corresponding to Q1 and Q3 profiles
        """

        #set up dimensions
        dimZ, dimX = np.shape(pm)     #get shape of data
        pmlvls = np.arange(0,dimZ*config.dz,config.dz)

        #locate centerline
        ctrZidx = pm.argmax(0)                          #locate maxima along height
        ctrXidx = pm.argmax(1)                          #locate maxima downwind
        pmCtr = np.array([pm[ctrZidx[nX],nX] for nX in range(dimX)])                #get concentration along the centerline

        xmax,ymax = np.nanargmax(ctrZidx), np.nanmax(ctrZidx)                       #get location of maximum centerline height
        centerline = ma.masked_where(pmlvls[ctrZidx] == 0, pmlvls[ctrZidx])         #make sure centerline is only calculated inside the plume
        centerline.mask[:int(1000/dx)] = True

        filter_window = max(int(utils.read_tag('W',[self.name])*10+1),51)
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

        #save attributes for quasi-stationary profile
        self.profile = interp1d(pmlvls,stableProfile,fill_value='extrapolate')(interpZ)
        self.quartiles = np.array([interpQ1,interpQ3])


        #calculate injection height variables ---------------------------
        zCL = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
        i_zCL = np.argmin(abs(interpZ - zCL))
        THzCL = self.sounding[i_zCL]

        self.zCL = zCL
        self.THzCL = THzCL
