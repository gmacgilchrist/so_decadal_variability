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

#################
# NOTE : nameposition is a super hacky solution to extract
# the name of the variable from the filename. It changes
# dependent on whether we are using the ocean grid or not.
# Need to overhaul procedure.

# Think I can combine these functions (to _get_specifics)
# but leaving separate for now
def _get_specifics_flux(fluxname):    
    specific={}
    specific['erai'] = {'suffix':'_1979-2018.nc',
                       'nameposition':slice(-24,-22)}
    specific['era5'] = {'suffix':'_1979-2019.nc',
                       'nameposition':slice(-20,-18)}
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
##################

def _get_specifics(name):    
    specific={}
    specific['erai'] = {'suffix':'_1979-2018.nc',
                       'nameposition':slice(-24,-22)}
    specific['era5'] = {'suffix':'_1979-2019.nc',
                       'nameposition':slice(-24,-22),
                       'gridlon':'longitude',
                       'gridlat':'latitude'}
    specific['jra55'] = {'suffix':'_1979-2019.nc',
                       'nameposition':slice(-25,-23)}
    specific['merra2'] = {'suffix':'_1980-2019.nc',
                       'nameposition':slice(-26,-24)}
    specific['en4'] = {'suffix':'_197901-201812.nc',
                      'depthname':'depth'}
    specific['iap'] = {'suffix':'_197901-201812.nc',
                      'depthname':'depth_std'}
    return specific[name]

## PATHS
def _get_oceanpath(oceanname, fluxname=None, varname=None):
    universal=_get_universal()
    specific=_get_specifics_ocean(oceanname)
    
    if fluxname is None:
        filename = universal['prefix']+'ocean_*'+oceanname+specific['suffix']
    else:
        filename = universal['prefix']+'ocean_*'+oceanname+'_'+fluxname+specific['suffix']
        
    if varname is not None:
        filename = universal['prefix']+'ocean_'+varname+'_'+oceanname+specific['suffix']
    path = universal['rootdir']+universal['localdir']+'ocean/'+filename
    return path

def _get_gridpath(name, varname=None):
    universal=_get_universal()
    specific=_get_specifics(name)
    
    filename = universal['prefix']+'grid_*'+name+universal['gridsuffix']
    
    if varname is not None:
        filename = universal['prefix']+'grid_'+varname+'_'+oceanname+universal['gridsuffix']
    path = universal['rootdir']+universal['localdir']+'grid/'+filename
    return path

def _get_fluxpath(fluxname, oceanname=None, varname=None):
    universal=_get_universal()
    specific=_get_specifics_flux(fluxname)
    if oceanname is None:
        filename = universal['prefix']+'flux_*'+fluxname+specific['suffix']
    else:
        filename = universal['prefix']+'flux_*'+fluxname+'_'+oceanname+specific['suffix']
    
    if varname is not None:
        filename = universal['prefix']+'flux_'+varname+'_'+fluxname+'_'+oceanname+specific['suffix']
    path = universal['rootdir']+universal['localdir']+'flux/'+filename
    return path

# LOADING
def _get_oceands(oceanname,fluxname=None):
    path = _get_oceanpath(oceanname,fluxname)
    return xr.open_mfdataset(path)

def _get_gridds(name):
    path = _get_gridpath(name)
    return xr.open_mfdataset(path)

def _get_fluxds(fluxname,oceanname=None):
    path = _get_fluxpath(fluxname,oceanname)
    return xr.open_mfdataset(path)
    
    # Getting flux data (more complicated loading because of muddled time coordinate)
    
    ### This was a hack that was in here before to select particular flux data
    ### I'm not entirely sure of its purpose, since the appropriate data seems to be
    ### collected irrespective. Will remove for now but may need to return to this.
    
#     universal=_get_universal()
#     specific=_get_specifics_flux(fluxname)
#     # Getting flux data (more complicated loading because of muddled time coordinate)
#     fluxfiles=[]
#     if oceanname is None:
#         filename = universal['prefix']+'flux_*'+fluxname+specific['suffix']
#     else:
#         filename = universal['prefix']+'flux_*'+fluxname+'_'+oceanname+specific['suffix']
#     for file in glob.glob(universal['rootdir']+'flux/'+filename):
#         print(file)
#         f = file[specific['nameposition']]
#         if f in ['sr','fw','ht']:
#             fluxfiles.append(file)
#     print('Retrieving data from :')
#     print(fluxfiles)
#     return xr.open_mfdataset(fluxfiles)

# PROCESSING
def _preprocess(fluxds,oceands,gridds,timeslice,onoceangrid):
    # HACK : current hack to avoid time selection for gridfile
    # when on flux grid (which is static), whereas ocean grid 
    # has time dimension
    # Revisit grid data to rectify
    timeselect = {'time':timeslice}
    if onoceangrid:
        fluxds = fluxds.sel(timeselect).assign_coords({'time':oceands['time'].sel(timeselect)})
        gridds = gridds.sel(timeselect)
        oceands = oceands.sel(timeselect)
    else:
        fluxds = fluxds.sel(timeselect)
#         gridds = gridds.sel(timeselect)
        oceands = oceands.sel(timeselect).assign_coords({'time':fluxds['time'].sel(timeselect)})
    # Merge
    ds = xr.merge([fluxds,oceands,gridds])
    # Roll longitude to it goes from 0 to 360
#     ds = ds.roll(lon=180,roll_coords=False).assign_coords({'lon':np.arange(0,360)})
    # Make heat flux positive into the ocean
    ds['ht'] *= -1
    ds['sr'] *= -1
    # Turn gamman to a proper density
    if 'gamman' in ds.data_vars:
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
def loaddata(fluxname, oceanname, timeslice, onoceangrid, debug=False):
    specific=_get_specifics_ocean(oceanname)
    if onoceangrid:
        # ocean
        oceands = _get_oceands(oceanname)
        # grid
        gridds = _get_gridds(oceanname)
        # flux
        fluxds = _get_fluxds(fluxname,oceanname)
    else:
        oceands = _get_oceands(oceanname,fluxname)
        gridds = _get_gridds(fluxname)
        fluxds = _get_fluxds(fluxname)
    
    if debug:
        return oceands,gridds,fluxds
    
    # Some renaming conventions
    if specific['depthname']!='depth':
        oceands = oceands.rename({specific['depthname']:'depth'})
        gridds = gridds.rename({specific['depthname']:'depth'})
        
    return _preprocess(fluxds,oceands,gridds,timeslice,onoceangrid)

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