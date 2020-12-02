import xarray as xr
import numpy as np
import glob

# PARAMETER SPECIFICATION
def _get_universal():
    universal={'rootdir': '/local/projects/so_decadal_variability/',
               'localdir': '',
               'prefix': 'SO_',
               'gridsuffix':'.nc'}
    return universal
    
def _get_specifics_flux(fluxname):    
    specific={}
    specific['erai'] = {'suffix':'_1979-2018.nc',
                       'nameposition':slice(-24,-22)}
    specific['era5'] = {'suffix':'_1979-2019.nc',
                       'nameposition':slice(-24,-22)}
    specific['jra55'] = {'suffix':'_1979-2019.nc',
                       'nameposition':slice(-25,-23)}
    specific['merra2'] = {'suffix':'_1980-2019.nc',
                       'nameposition':slice(-26,-24)}
    return specific[fluxname]

def _get_specifics_ocean(oceanname):
    specific={}
    specific['en4'] = {'suffix':'_197901-201812.nc',
                      'depthname':'depth'}
    specific['iap'] = {'suffix':'_197901-201812.nc',
                      'depthname':'depth_std'}
    return specific[oceanname]

## PATHS
def _get_oceanpath(oceanname, varname=None):
    universal=_get_universal()
    specific=_get_specifics_ocean(oceanname)
    
    filename = universal['prefix']+'ocean_*'+oceanname+specific['suffix']
    if varname is not None:
        filename = universal['prefix']+'ocean_'+varname+'_'+oceanname+specific['suffix']
    path = universal['rootdir']+universal['localdir']+'ocean/'+filename
    return path

def _get_gridpath(oceanname, varname=None):
    universal=_get_universal()
    specific=_get_specifics_ocean(oceanname)
    
    filename = universal['prefix']+'grid_*'+oceanname+universal['gridsuffix']
    if varname is not None:
        filename = universal['prefix']+'grid_'+varname+'_'+oceanname+universal['gridsuffix']
    path = universal['rootdir']+universal['localdir']+'grid/'+filename
    return path

def _get_fluxpath(oceanname, fluxname, varname=None):
    universal=_get_universal()
    specific=_get_specifics_flux(fluxname)
    filename = universal['prefix']+'flux_*'+fluxname+'_'+oceanname+specific['suffix']
    if varname is not None:
        filename = universal['prefix']+'flux_'+varname+'_'+fluxname+'_'+oceanname+specific['suffix']
    path = universal['rootdir']+universal['localdir']+'flux/'+filename
    return path

# LOADING
def _get_oceands(oceanname):
    path = _get_oceanpath(oceanname)
    return xr.open_mfdataset(path)

def _get_gridds(oceanname):
    path = _get_gridpath(oceanname)
    return xr.open_mfdataset(path)

def _get_fluxds(fluxname,oceanname):
    universal=_get_universal()
    specific=_get_specifics_flux(fluxname)
    # Getting flux data (more complicated loading because of muddled time coordinate)
    fluxfiles=[]
    filename = universal['prefix']+'flux_*'+fluxname+'_'+oceanname+specific['suffix']
    for file in glob.glob(universal['rootdir']+'flux/'+filename):
        f = file[specific['nameposition']]
        if f in ['sr','fw','ht']:
            fluxfiles.append(file)
    return xr.open_mfdataset(fluxfiles)

# PROCESSING
def _preprocess(fluxds,oceands,gridds,timeslice):
    timeselect = {'time':timeslice}
    fluxds = fluxds.sel(timeselect).assign_coords({'time':oceands['time'].sel(timeselect)})
    gridds = gridds.sel(timeselect)
    oceands = oceands.sel(timeselect)
    # Merge
    ds = xr.merge([fluxds,oceands,gridds])
    # Roll longitude to it goes from 0 to 360
#     ds = ds.roll(lon=180,roll_coords=False).assign_coords({'lon':np.arange(0,360)})
    # Make heat flux positive into the ocean
    ds['ht'] *= -1
    ds['sr'] *= -1
    # Turn gamman to a proper density
    ds['gamman']+=1000
    return ds

def _preprocess_oceanonly(oceands,gridds,timeslice, roll):
    timeselect = {'time':timeslice}
    gridds = gridds.sel(timeselect)
    oceands = oceands.sel(timeselect)
    # Merge
    ds = xr.merge([oceands,gridds])
    # Roll longitude to it goes from 0 to 360
    if roll:
        ds = ds.roll(lon=180,roll_coords=False).assign_coords({'lon':np.arange(0,360)})
    # Turn gamman to a proper density
    ds['gamman']+=1000
    return ds

# LOADING WRAPPERS
def loaddata(fluxname, oceanname, timeslice, debug=False):
    specific=_get_specifics_ocean(oceanname)
    # ocean
    oceands = _get_oceands(oceanname)
    # grid
    gridds = _get_gridds(oceanname)
    # flux
    fluxds = _get_fluxds(fluxname,oceanname)
    
    if debug:
        return oceands,gridds,fluxds
    
    # Some renaming conventions
    if specific['depthname']!='depth':
        oceands = oceands.rename({specific['depthname']:'depth'})
        gridds = gridds.rename({specific['depthname']:'depth'})
        
    return _preprocess(fluxds,oceands,gridds,timeslice)

def loaddata_oceanonly(oceanname, timeslice, roll=True):
    
    specific=_get_specifics_ocean(oceanname)
    
    oceands = _get_oceands(oceanname)
    gridds = _get_gridds(oceanname)
    
    # Some renaming conventions
    if specific['depthname']!='depth':
        oceands = oceands.rename({specific['depthname']:'depth'})
        gridds = gridds.rename({specific['depthname']:'depth'})
        
    return _preprocess_oceanonly(oceands,gridds,timeslice, roll)

# POST-PROCESSING AND SAVING
def save_ocean(da,oceanname):
    universal = _get_universal()
    specific = _get_specifics_ocean(oceanname)
    # Change naming conventions back
    if specific['depthname']!='depth':
        da = da.rename({'depth':specific['depthname']})
    varname = da.name
    filename = universal['prefix']+'ocean_'+varname+'_'+oceanname+specific['suffix']
    path = universal['rootdir']+universal['localdir']+'ocean/'+filename
    print('Saving to '+path)
    da.to_netcdf(path)