import xarray as xr
import numpy as np
from xhistogram.xarray import histogram

# Watermass transformation calculation
def _calc_shortwave_penetration(ds,xgrid):
    R = 0.58
    h1 = 0.35
    h2 = 23

    # Calculate in 1D
    Fps77 = R*np.exp(-ds['depth_i']/h1)+(1-R)*np.exp(-ds['depth_i']/h2)

    # Now get a 4D shortwave field
    shortwave = ds['sr']*Fps77
    # and take the divergence
    # (note that for taking the derivative, 'depth' must be positive upward, so reverse the sign)
    # (subsequently, reverse the sign so that it is consistent with the other variables)
    dsr4d = xgrid.derivative(-shortwave,'Z')
    
    return dsr4d

def _calc_densityflux(FW,Q,Qsw,S,alpha,beta,Cp=4200):
    
    F = xr.Dataset()
    F['heat'] = (alpha/Cp)*Q + (alpha/Cp)*Qsw
    F['fw'] = -FW*S*beta
    F['total'] = F['heat']+F['fw']
    
    return F

def calc_densityflux(ds,xgrid):
    dsr4d = _calc_shortwave_penetration(ds,xgrid)
    mask = _create_mask(ds)
    F = _calc_densityflux(ds['fw']*mask,ds['ht']*mask,dsr4d,ds['sa'],ds['alpha'],ds['beta'])
    
    return F

def _create_mask(ds):
    # Create a 3D mask with 1/dz in the surface and zero elsewhere
    idz = 1/ds['dz'][0].values
    mask = xr.concat([idz*xr.ones_like(ds['ct'].isel(time=0,depth=0)),
                      xr.zeros_like(ds['ct'].isel(time=0,depth=slice(1,None)))],
                     dim='depth')
    return mask

### Watermass transformation calculation
    
def _calc_watermasstransformation(F,gamman,b,V,gn_edges):
    # Discrete volume calculation derived in Appendix 7.5 of Groeskamp et al (2018)
    G = xr.Dataset()
    for var in F.data_vars:
        gFbV = gamman*b*F[var]*V
        nanmask=np.isnan(gFbV)
        G[var] = histogram(gamman.where(~nanmask),bins=[gn_edges],weights=gFbV.where(~nanmask),dim=['lat','lon','depth'])/np.diff(gn_edges)
    return G

def calc_watermasstransformation(ds,xgrid,gn_edges):
    dsr4d = _calc_shortwave_penetration(ds,xgrid)
    mask = _create_mask(ds)
    F = _calc_densityflux(ds['fw']*mask,ds['ht']*mask,dsr4d,ds['sa'],ds['alpha'],ds['beta'],Cp=4200)
    G = _calc_watermasstransformation(F,ds['gamman'],ds['b'],ds['vol4d'],gn_edges)
    
    return G

### Storage change calculation

def _calc_dMdt(mass,gamman,gn_edges):
    # Augment edges
    gn_edges_all = np.concatenate(([np.array(-99999)],gn_edges,[np.array(99999)]))
    # Histogram mass
    nanmask=np.isnan(gamman)
    M_on_gamma = histogram(gamman.where(~nanmask),
                       bins=[gn_edges_all],
                       weights=mass.where(~nanmask),
                       dim=['lat','lon','depth']).transpose()
    
    # To integrate for all volume with temperature greater than a certain value,
    # take cumulative sum and reassign the coordinates to align with G
    M_on_gamma_cumsum = xr.concat([xr.zeros_like(M_on_gamma.isel({'gamman_bin':0})),
                                   M_on_gamma.cumsum('gamman_bin')],dim='gamman_bin')
    # We wish to have the total mass for the volume with temperature greater than that contour,
    # So take away the total sum from the cumulative sum to reverse the direction
    M_reverse = (M_on_gamma.sum('gamman_bin')-M_on_gamma_cumsum)
    # Now we can get rid of the boundary contours, which were there to ensure that all
    # of the volume wass captures, and we assign the coordinates to match with G
    M = M_reverse.isel(gamman_bin=slice(1,-1)).assign_coords({'gamman_bin':gn_edges})

    # Calculate the derivative with respect to time
    dMdt = M.diff('time')/(M['time'].diff('time').astype('float')*1E-9)
    
    # The time derivative is align with the start of each month,
    # so define a new time coordinate
    timenew = M.time[:-1]+(M['time'].shift({'time':-1})-M['time'][:-1])/2
    
    # Assign that coordinate for the time derivative
    dMdt = dMdt.assign_coords(time=timenew)

    # Rename
    dMdt.name = 'dMdt'
    
    return dMdt

def calc_dMdt(ds,gn_edges):
    return _calc_dMdt(ds['mass'],ds['gamman'],gn_edges)

### b-factor

# Wrappers for derivative operations in xgcm
def _xgcm_interp_and_derivative(da,xgrid,dim,boundary=None):
    # Interpolate to grid cell boundaries
    da_i = xgrid.interp(da,dim,boundary=boundary)
    # Take the derivative
    dadl = xgrid.derivative(da_i,dim,boundary=boundary)
    return dadl
    
def _xgcm_interp_and_derivative_3D(da,xgrid,dims=['X','Y','Z'],boundaries=[None,None,None]):
    
    # Calculate gradients in X, Y and Z
    dad1 = _xgcm_interp_and_derivative(da,xgrid,dims[0],boundaries[0])
    dad2 = _xgcm_interp_and_derivative(da,xgrid,dims[1],boundaries[1])
    dad3 = _xgcm_interp_and_derivative(da,xgrid,dims[2],boundaries[2])
    
    return dad1, dad2, dad3

def _calc_bfactor(T,S,rho,alpha,beta,gamma,xgrid):
    
    # Derivatves in T, S, and gamma
    dims = ['X','Y','Z']
    boundaries = [None,'extend','extend']
    dTdx,dTdy,dTdz = _xgcm_interp_and_derivative_3D(T,xgrid,dims,boundaries)
    dSdx,dSdy,dSdz = _xgcm_interp_and_derivative_3D(S,xgrid,dims,boundaries)
    dgdx,dgdy,dgdz = _xgcm_interp_and_derivative_3D(gamma,xgrid,dims,boundaries)
    
    # Locally referenced potential density
    drdx = rho*(-alpha*dTdx + beta*dSdx)
    drdy = rho*(-alpha*dTdy + beta*dSdy)
    drdz = rho*(-alpha*dTdz + beta*dSdz)

    # Calculate the absolute magnitudes
    abs_drd = xr.ufuncs.sqrt(xr.ufuncs.square(drdx)+xr.ufuncs.square(drdy)+xr.ufuncs.square(drdz))
    abs_dgd = xr.ufuncs.sqrt(xr.ufuncs.square(dgdx)+xr.ufuncs.square(dgdy)+xr.ufuncs.square(dgdz))
        
    # Calculate ratio
    b = abs_drd/abs_dgd
    
    return b

def calc_bfactor(ds,xgrid):
    return _calc_bfactor(ds['ct'],ds['sa'],ds['rho'],ds['alpha'],ds['beta'],ds['gamman'],xgrid)