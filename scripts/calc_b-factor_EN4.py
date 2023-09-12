# Python script to calculate the b-factor from EN4 hydrographic data
# gmac 8/10/20
# gmac 2/12/20: adjusting to also work with iap

import xarray as xr
import numpy as np
from xgcm import Grid
from numba import jit
import glob

rootdir = '/local/projects/so_decadal_variability/'
localdir = ''
prefix = 'SO_'
grid_name = 'iap'
ocean_name = 'iap'
grid_suffix = '.nc'
ocean_suffix = '_197901-201812.nc'

# Ocean data
filename = prefix+'ocean_*'+ocean_name+ocean_suffix
ocean = xr.open_mfdataset(rootdir+localdir+'ocean/'+filename)
# Grid data
filename = prefix+'grid_*'+grid_name+grid_suffix
grid = xr.open_mfdataset(rootdir+'grid/*'+filename)
ds = xr.merge([ocean,grid])

if ocean_name=='iap':
    ds = ds.rename({'depth_std':'depth'})

# Create xgcm Grid object
@jit(nopython=True)
def haversine_distance(lat1,lat2,lon1,lon2,degrees=True):
    R=6371E3
    if degrees:
        fac = np.pi/180
        lat1 = lat1*fac
        lat2 = lat2*fac
        lon1 = lon1*fac
        lon2 = lon2*fac
    return 2*R*np.arcsin(np.sqrt(np.sin((lat2-lat1)/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2))

@jit(nopython=True)
def calc_spacing_from_latlon(lat,lon):
    ni = len(lon)-1
    nj = len(lat)-1
    dx = np.empty(shape=(ni,nj))
    dy = np.empty(shape=(ni,nj))
    for i in range(ni):
        lon1 = lon[i]
        lon2 = lon[i+1]
        for j in range(nj):
            lat1 = lat[j]
            lat2 = lat[j+1]

            dx[i,j] = haversine_distance(lat1,
                                         lat1,
                                         lon1,
                                         lon2,
                                         degrees=True)
            dy[i,j] = haversine_distance(lat1,
                                         lat2,
                                         lon1,
                                         lon1,
                                         degrees=True)
    return dx, dy

def _calc_outerdepths(depth):
    '''Given the cell-centre depth points, determine the outer points.
    Assumption that the first outer point is at the ocean surface, z=0.'''
    
    nk = len(depth)
    depth_i_vals = np.zeros(nk+1)
    for k in range(nk):
        if k>0:
            depth_i_vals[k] = (depth[k-1]+depth[k])/2
    depth_i_vals[nk] = depth[nk-1]+(depth[nk-1]-depth_i_vals[nk-1])
    return depth_i_vals

# Specify coordinates with ghost points
lonG = np.append(ds['lon']+0.5,ds['lon'][-1]+1.5,)
latG = np.append(ds['lat']+0.5,ds['lat'][-1]+1.5,)
latC = np.append(ds['lat'],ds['lat'][-1]+1)
lonC = np.append(ds['lon'],ds['lon'][-1]+1)
# Find distances
# ... across cell faces
dxC,_ = calc_spacing_from_latlon(latC,lonG) # defined at latC,lonC
_,dyC = calc_spacing_from_latlon(latG,lonC) # defined at latC,lonC
# ... between grid cell boundaries
dxG,dyG = calc_spacing_from_latlon(latC,lonC) # defined at latC,lonG & latG,lonC
# Specify cell boundary coordinates
ds['lonG'] = xr.DataArray(lonG[:-1],dims=['lonG'])
ds['latG'] = xr.DataArray(latG[:-1],dims=['latG'])
# Specify grid spacings
ds['dxC'] = xr.DataArray(dxC,dims=['lon','lat'])
ds['dyC'] = xr.DataArray(dyC,dims=['lon','lat'])
ds['dxG'] = xr.DataArray(dxG,dims=['lonG','lat'])
ds['dyG'] = xr.DataArray(dyG,dims=['lon','latG'])
# Find depths of outer points
if 'depth_bnds' in ds.data_vars:
    depth_i_vals = np.append(grid.depth_bnds.isel(bnds=0), grid.depth_bnds.isel(bnds=1)[-1])
else:
    depth_i_vals = _calc_outerdepths(ds['depth'])
ds['depth_i'] = xr.DataArray(depth_i_vals,coords={'depth_i':depth_i_vals},dims={'depth_i'})

# Get depth distance
ds['dz'] = ds['depth_i'].diff('depth_i').rename({'depth_i':'depth'}).assign_coords({'depth':ds['depth']})

coords = {
    'X':{'center':'lon','right':'lonG'},
    'Y':{'center':'lat','right':'latG'},
    'Z':{'center':'depth','outer':'depth_i'}
}
metrics = {
    'X':['dxC','dxG'],
    'Y':['dyC','dyG'],
    'Z':['dz']
}
xgrid = Grid(ds,coords=coords,metrics=metrics,periodic=['X'])

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

def calc_b(T,S,rho,alpha,beta,gamma,xgrid):
    
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

# Calculate b
b = calc_b(ds['ct'],ds['sa'],ds['rho'],ds['alpha'],ds['beta'],ds['gamman'],xgrid)

# Save
print('Saving netcdf...')
filename = prefix+'ocean_b_'+ocean_name+ocean_suffix
b.name = 'b'
b.to_netcdf(rootdir+localdir+'ocean/'+filename)