# Define an xgcm grid object given one-dimensional lat, lon, depth
# coordinates for the EN4 dataset

# gmac 8/10/20

# Note that this is not actually an executable script at the moment
# I am simply placeholding the functions and flow, with a view to 
# creating a small wrapper package to perform a load of these functions

import xarray as xr
import numpy as np
from xgcm import Grid
from numba import jit

# Load data
rootdir = '/local/projects/so_decadal_variability/'
localdir = ''
prefix = 'SO_'
grid_name = 'en4'
grid_suffix = '.nc'
# Grid data
filename = prefix+'grid_*'+grid_name+grid_suffix
grid = xr.open_mfdataset(rootdir+'grid/*'+filename)
ds = xr.merge([ocean,flux,grid])

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
depth_i_vals = np.append(grid.depth_bnds.isel(bnds=0), grid.depth_bnds.isel(bnds=1)[-1])
ds['depth_i'] = xr.DataArray(depth_i_vals,coords={'depth_i':depth_i_vals},dims={'depth_i'})
# Get depth distance
ds['dz'] = ds['depth_i'].diff('depth_i').rename({'depth_i':'depth'}).assign_coords({'depth':ds['depth']})

# Define xgcm object
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