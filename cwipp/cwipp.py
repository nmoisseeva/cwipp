import config
import utils
import graphics
import numpy as np
from numpy import ma
from scipy.optimize import fsolve
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


class Plume:
    """
    Parent Plume class.

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
    wf : float
        characteristic fire velocity scale [m s-1]
    Tau : float
        characteristic timescale [s]

    """

    def __init__(self, name):
        """
        Constructs the plume object with some inital attributes

        Parameters
        -----------
        name: str
            plume name
        """
        #
        # #get initial raw sounding (from cross-section wrfcs data)
        # T0 = np.load(config.wrfdir + 'profiles/profT0' + name + '.npy')    #load initial temperature profile
        #
        # #get BL height
        # zi = utils.get_zi(T0,config.dz)                    #calculate BL height
        # zs = zi * config.BLfrac
        #
        # #interpolate sounding to analysis levels
        # metlvls = np.arange(0,len(T0)*config.dz,config.dz)
        # interpT= interp1d(metlvls,T0,fill_value='extrapolate')
        # T0interp = interpT(config.interpZ)
        # i_zs = np.argmin(abs(config.interpZ - zs))
        # THs = T0interp[i_zs]

        self.name = name
        # self.zi = zi
        # self.zs = zs
        # self.sounding = T0interp
        # self.THs = THs

    def get_sounding(self, T0):
        """
        Calculates attributes relating to vertical potential temperature profile

        Parameters
        -----------
        T0: ndarray
            potential temperature profile on host model levels (not interpolated)

        Returns
        ---------
        zi : float
            boundary layer height [m]
        zs : float
            refernce height (zi * BLfrac) [m]
        sounding: ndarray
            vertical potential temperature sounding on interpolated analysis levels [K]
        THs : float
            ambient potential tempreature at reference height zs [K]
        """

        # #get initial raw sounding (from cross-section wrfcs data)
        # T0 = np.load(config.wrfdir + 'profiles/profT0' + name + '.npy')    #load initial temperature profile

        #get BL height
        zi = utils.get_zi(T0,config.dz)                    #calculate BL height
        zs = zi * config.BLfrac

        #interpolate sounding to analysis levels
        metlvls = np.arange(0,len(T0)*config.dz,config.dz)
        interpT= interp1d(metlvls,T0,fill_value='extrapolate')
        T0interp = interpT(config.interpZ)
        i_zs = np.argmin(abs(config.interpZ - zs))
        THs = T0interp[i_zs]

        self.zi = zi
        self.zs = zs
        self.sounding = T0interp
        self.THs = THs


    def get_wf(self):
        """

        Finds characteristic time (Tau) and velocity (wf) scales.

        Returns
        ---------
        wf : float
            characteristic fire velocity scale [m s-1]
        Tau : float
            characteristic timescale [s]
        """
        Tau = 1/np.sqrt(config.g*(self.THzCL - self.THs)/(self.THs * (self.zCL-self.zs)))
        wf= ((config.g*self.I*(self.zCL-self.zs))/(self.THs*self.zi))**(1/3.)

        self.Tau = Tau
        self.wf = wf


    def classify(self):
        """

        Classifies the plume as penetrative (True) or boundary layer (False)

        Returns
        --------
        penetrative : boolean
            classification (True if penetrative).
        """

        if  self.zCL < (self.zi + (config.dz)/2):
            self.penetrative = False
        else:
            self.penetrative = True


class LESplume(Plume):
    """
    Child Plume class used for simulated plumes (non-predictive mode).

    Assumes full model fields are available.

    Attributes
    ----------
    profile : ndarray
        1D vector corresponding to quasi-stationary downwind PM profile [concentration]
    quartiles : ndarray
        2D array with columns corresponding to Q1 and Q3 profiles [concentration]
    I : float
        fireline intensity parameter [K m2 s-1]
    zCL : float
        plume injection height [m]
    centerline: ndarray
        masked array containing height indices of plume centerline
    ctr_idx: list
        list of vertical indecies correponding to centerline height
    THzCL : float
        ambient potential temperature at zCL [K]

    """

    def get_I(self, flux, length, *Us):
        """
        Finds cross-wind fireline intensity parameter I

        Parameters
        -----------
        flux : ndarray
            3D (time,y,x) array containing heat flux values [kW m-2].
        length : float
            maximum cross-wind length of the fire over the entire timespan [m].
        Us : float, optional
            surface wind direction [deg, relative to y axis] NOT CURRENTLY IMPLEMENTED!

        Returns
        ---------
        I : float
            fireline intensity parameter [K m2 s-1]

        """

        #confirm there are sufficiant dimensions
        dims = np.shape(flux)
        if len(dims) > 3:
            raise ValueError('Too many dimensions for heat flux data')
        elif len(dims)<3:
            raise ValueError('Too few dimensions: need 3D array (time,y,x)')

        #mask and pad the heat source ------------------------
        upwind_padding = int(length/config.dx)
        downwind_padding = int(2000/config.dx)              #assumes ground is not heated beyont 1km downwind
        masked_flux = ma.masked_less_equal(np.pad(flux,((0,0),(0,0),(upwind_padding,0)), 'constant',constant_values=0),1)

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

    def get_zCL(self, pm, **kwargs):
        """
        Extracts mean injection height from LES.

        Finds quasi-stationary downwind profile and its IQR, extracts injection height and associated variables.

        Parameters
        ----------
        pm : ndarray
            2D array (z,x) of pm cross-section
        plot: boolean, optional
            create a multi-panel plot of the method, requires csdict argument to follow
        csdict: dict, optional
            cross-section dictionary for the plume, if plotting is required

        Returns
        --------
        profile : ndarray
            1D vector corresponding to quasi-stationary downwind PM profile
        quartiles: ndarray
            2D array with columns corresponding to Q1 and Q3 profiles
        zCL : float
            plume injection height [m]
        centerline: ndarray
            masked array containing height indices of plume centerline
        ctr_idx: list
            list of vertical indecies correponding to centerline height
        THzCL : float
            ambient potential temperature at zCL [K]
        """

        import warnings
        warnings.filterwarnings("ignore")

        #set up dimensions
        dimZ, dimX = np.shape(pm)     #get shape of data
        pmlvls = np.arange(0,dimZ*config.dz,config.dz)

        #locate centerline
        ctrZidx = np.nanargmax(pm,0)                          #locate maxima along height
        ctrXidx = np.nanargmax(pm,1)                          #locate maxima downwind
        i_zi = np.nanargmin(abs(pmlvls - self.zi))

        ctr_idx = []
        for nX in range(dimX):
            if nX <  ctrXidx[0]:
                idx = 0
            # elif nX < ctrXidx[i_zi]:
            elif ctr_idx[-1]<i_zi-1:
                idx = np.nanargmax(pm[:i_zi,:],0)[nX]
                closestZ = np.nanargmin(abs(ctrXidx - nX))
                if idx > closestZ or idx==0:
                    if closestZ < i_zi:
                        idx = closestZ
            else:
                idx = ctrZidx[nX]
            ctr_idx.append(idx)


        PMctr = np.array([pm[ctr_idx[nX],nX] for nX in range(dimX)])                #get concentration along the centerline

        xmax,ymax = np.nanargmax(ctr_idx), np.nanmax(ctr_idx)                       #get location of maximum centerline height
        centerline = ma.masked_where(pmlvls[ctr_idx] == 0, pmlvls[ctr_idx])         #make sure centerline is only calculated inside the plume
        centerline.mask[:int(1000/config.dx)] = True

        filter_window = max(int(utils.read_tag('W',[self.name])*10+1),51)
        smoothCenterline = savgol_filter(centerline, filter_window, 3)              #smooth centerline height

        #calculate concentration changes along the centerline
        dPMdX = PMctr[1:]-PMctr[0:-1]
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

        if sum(stablePMmask) == 0:
            stablePMmask = [True if abs(smoothPM[nX])< np.nanmax(smoothPM)*0.1 and \
                                    nX > np.nanargmax(centerline[~centerline.mask][:-50]) and\
                                    nX > np.nanargmax(smoothPM) else\
                                    False for nX in range(dimX-1) ]
        stablePM = pm[:,1:][:,stablePMmask]
        stableProfile = np.median(stablePM,1)
        #find IQR
        pmQ1 = np.percentile(stablePM,25,axis = 1)
        pmQ3 = np.percentile(stablePM,75,axis = 1)
        interpQ1 = interp1d(pmlvls,pmQ1,fill_value='extrapolate')(config.interpZ)
        interpQ3 = interp1d(pmlvls,pmQ3,fill_value='extrapolate')(config.interpZ)

        #save attributes for quasi-stationary profile
        self.profile = interp1d(pmlvls,stableProfile,fill_value='extrapolate')(config.interpZ)
        self.quartiles = np.array([interpQ1,interpQ3])
        self.centerline = centerline
        self.ctr_idx = ctr_idx

        #calculate injection height variables
        zCL = np.mean(smoothCenterline[1:][stablePMmask])    #injection height is where the centerline is stable and concentration doesn't change
        i_zCL = np.argmin(abs(config.interpZ - zCL))
        THzCL = self.sounding[i_zCL]

        self.zCL = zCL
        self.THzCL = THzCL
        #
        #make plots, if requested
        if len(kwargs.keys()) > 0:
            if kwargs['plot']:
                fireCS = kwargs['csdict']['ghfx'][-1,:]
                flux2D = kwargs['csdict']['ghfx2D'][-1,:,:]
                graphics.plot_zcl(self,pm,fireCS,flux2D,stablePMmask,smoothCenterline)

class MODplume(Plume):
    """
    Child Plume class used for modelled plumes (predictive mode)

    Attributes
    ----------
    zCL : float
        parameterized plume injection height [m]
    THzCL : float
        ambient potential temperature at modelled zCL [K]

    Methods
    -------
    iterate(self):
        Applies iterative solution to parameterize plume injection height
    """

    def iterate(self, biasFit=None, **kwargs):
        """
        Applies iterative solution to parameterize plume injection height

        Parameters
        ----------
        biasFit : array_like, optional
            bias fit parameters. If none provided defaults to m = 1, b = 0.
        argout: boolean, optional
            flag to output return arguments. If False(default) they are assigned as attributes

        Returns
        -------
        zCL : float
            parameterized plume injection height [m]
        THzCL : float
            ambient potential temperature at modelled zCL [K]
        """
        if biasFit:
            m, b = biasFit[0], biasFit[1]
        else:
            m, b = 1, 0

        i_zs = np.nanargmin(abs(config.interpZ - self.zs))

        toSolve = lambda z : z  - b - m*(self.zs + \
                        1/(np.sqrt(config.g*(self.sounding[int(z/config.zstep)] - self.THs)/(self.THs * (z-self.zs))))  * \
                        (config.g*self.I*(z-self.zs)/(self.THs * self.zi))**(1/3.))

        zCL = fsolve(toSolve, self.zi, factor=0.1)
        i_zCL = np.nanargmin(abs(config.interpZ - zCL))

        THzCL = self.sounding[i_zCL]

        if 'argout' in kwargs.keys():
            if kwargs['argout']:
                return float(zCL), THzCL
            else:
                self.THzCL = THzCL
                self.zCL = float(zCL)
        else:
            self.THzCL = THzCL
            self.zCL = float(zCL)

    def explicit_solution(self, Gamma, ze, biasFit=None):
        """
        Applies explicit solution to parameterize plume injection height

        Parameters
        ----------
        biasFit : array_like, optional
            bias fit parameters. Default is m = 1, b = 0

        Returns
        -------
        zCL : float
            parameterized plume injection height [m]
        THzCL : float
            ambient potential temperature at modelled zCL [K]
        """
        if biasFit:
            m, b = biasFit[0], biasFit[1]
        else:
            m, b = 1, 0

        zCL = m*(((self.THs/config.g)**(1/4.)) * ((self.I/self.zi)**(0.5)) * ((1/Gamma)**(3/4.)) + ze) + b

        self.zCL = zCL

'''
FUTURE DEVELOPMENT:
    def get_profile(self):
        """
        Parameterization of the full normalized vertical smoke profile


        Parameters
        ----------


        Returns
        -------
        profile : ndarray
            1D vector corresponding to quasi-stationary downwind PM profile
        """

        profile = np.empty((len(config.interpZ))) * np.nan

        if not self.penetrative:
            profile = 1
        elif self.penetrative:
            self.get_wf()

            #get Deadorff's velocity for spread
            wD = (config.g * self.zi * 0.13 / self.THs)**(1/3.) #!!!! HARDCODED SURFACE HEAT FLUX

            sigma_top = (self.zCL - self.zs)/3.
            if self.wf/wD < 1.5:
                Rw = self.U/self.wf
            else:
                Rw = self.U/(self.wf - wD)

            if Rw > 1:
                sigma_bottom = Rw * sigma_top
            else:
                sigma_bottom = sigma_top

            izCL = np.argmin(abs(config.interpZ - self.zCL))
            profile[izCL:] = np.exp(-0.5*((config.interpZ[izCL:] - self.zCL)/sigma_top)**2)
            profile[:izCL+1] = np.exp(-0.5*((config.interpZ[:izCL+1] - self.zCL)/sigma_bottom)**2)

            self.profile = profile
'''
