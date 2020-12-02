import xarray as xr
import numpy as np
import glob

# Loading and processing data
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

def _get_oceands(oceanname):
    universal=_get_universal()
    specific=_get_specifics_ocean(oceanname)
    filename = universal['prefix']+'ocean_*'+oceanname+specific['suffix']
    return xr.open_mfdataset(universal['rootdir']
                              +universal['localdir']
                              +'ocean/'
                              +filename)
def _get_gridds(oceanname):
    universal=_get_universal()
    specific=_get_specifics_ocean(oceanname)
    filename = universal['prefix']+'grid_*'+oceanname+universal['gridsuffix']
    return xr.open_mfdataset(universal['rootdir']
                              +universal['localdir']
                              +'grid/'
                              +filename)

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

def _preprocess(fluxds,oceands,gridds,timeslice):
    timeselect = {'time':timeslice}
    fluxds = fluxds.sel(timeselect).assign_coords({'time':oceands['time'].sel(timeselect)})
    gridds = gridds.sel(timeselect)
    oceands = oceands.sel(timeselect)
    # Merge
    ds = xr.merge([fluxds,oceands,gridds])
    # Roll longitude to it goes from 0 to 360
    ds = ds.roll(lon=180,roll_coords=False).assign_coords({'lon':np.arange(0,360)})
    # Make heat flux positive into the ocean
    ds['ht'] *= -1
    ds['sr'] *= -1
    # Turn gamman to a proper density
    ds['gamman']+=1000
    return ds

def _preprocess_oceanonly(oceands,gridds,timeslice):
    timeselect = {'time':timeslice}
    gridds = gridds.sel(timeselect)
    oceands = oceands.sel(timeselect)
    # Merge
    ds = xr.merge([oceands,gridds])
    # Roll longitude to it goes from 0 to 360
    ds = ds.roll(lon=180,roll_coords=False).assign_coords({'lon':np.arange(0,360)})
    # Turn gamman to a proper density
    ds['gamman']+=1000
    return ds

def loaddata(fluxname, oceanname, timeslice):
    specific=_get_specifics_ocean(oceanname)
    # ocean
    oceands = _get_oceands(oceanname)
    # grid
    gridds = _get_gridds(oceanname)
    # flux
    fluxds = _get_fluxds(fluxname,oceanname)
    
    # Some renaming conventions
    if specific['depthname']!='depth':
        oceands = oceands.rename({specific['depthname']:'depth'})
        gridds = gridds.rename({specific['depthname']:'depth'})
        
    return _preprocess(fluxds,oceands,gridds,timeslice)

def loaddata_oceanonly(oceanname, timeslice):
    oceands = _get_oceands(oceanname)
    gridds = _get_gridds(oceanname)
    return _preprocess_oceanonly(oceands,gridds,timeslice)