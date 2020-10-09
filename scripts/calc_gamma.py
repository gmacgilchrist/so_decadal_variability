#!/nbhome/gam/miniconda/envs/mom6-clean/bin/python

# Calculate gamma_n from Jackett and McDougall (1997) using 
# the python wrapper developed by E. Firing 
# (https://currents.soest.hawaii.edu/hgstage/pygamma/)
print('---        ALLEZ        ---')

print('---  Importing modules  ---')
import xarray as xr
from pygamma import gamma_n
import wmt_bgc.basic as wmt
import numpy as np
import glob
import re
from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return dict(islice(iterable, n))

# Attributes for gamma_n
g_attrs = {'long_name':'Neutral density (gamma_n) from Jackett and McDougall (1997)',
           'units':'kgm-3',
           'standard_name':'gamma_n'}

# gamma_n function for 4-D data
def _gamma_n_4D(S,T,p,lon,lat,verbose=False):
    ni,nj,nk,nt = S.shape
    g_vals = np.empty(shape=(ni,nj,nk,nt))
    for t in range(nt):
        for j in range(nj):
            if verbose:
                if np.mod(j,100)==0:
                    print('time-index = '+str(t)+'; y-index = '+str(j))
            g_vals[:,j,:,t],_,_ = gamma_n(S[:,j,:,t],
                                          T[:,j,:,t],
                                          p[:,j,:,t],
                                          lon[:,j],
                                          lat[:,j])
    return g_vals

# Define paths
rootdir = '/archive/Raphael.Dussin/xanadu_esm4_20190304_mom6_2019.08.08/OM4p25_JRA55do1.4_0netfw_cycle6/gfdl.ncrc4-intel16-prod/pp/'
pp = 'ocean_annual_z'
localdir = '/ts/annual/5yr/'
filename = pp+'.*.'
files = glob.glob(rootdir+pp+localdir+filename+'so.nc')
# Load grid
grid = xr.open_dataset(rootdir+pp+'/'+pp+'.static.nc')[['geolon','geolat']]

outdir = '/archive/gam/so_decadal_variability/OM4p25_JRA55do1.4_0netfw_cycle6/'

# Reconfigure longitude to be 0 to 360
grid['geolon'] = grid['geolon'].where(grid['geolon']>0,grid['geolon']+360)
# Load depth levels
grid['z_l']=xr.open_dataset(files[0])['z_l']
# Calculate pressure
grid['p'] = wmt.gsw_p_from_z(-1*grid['z_l'],grid['geolat'])

print('--- Looping through years ---')
for file in files:
    # Get the year range
    m = re.search(localdir+pp+'.(.+?).so.nc', file)
    print(m.group(1))
    inpath = rootdir+pp+localdir+pp+'.'+m.group(1)+'.'
    outpath = outdir+pp+'/'+pp+'.'+m.group(1)+'.'
    
    ds = xr.Dataset()
    ds['thetao'] = xr.open_dataset(inpath+'thetao.nc')['thetao']
    ds['so'] = xr.open_dataset(inpath+'so.nc')['so']
    
    select = {'xh':slice(500,510),'yh':slice(200,210),'z_l':slice(0,5),'time':slice(0,2)}
    dsnow = ds.transpose()#.isel(select)
    gridnow = grid.transpose()#.isel(take(3,select.items()))

    g_vals = _gamma_n_4D(dsnow['so'],
                         dsnow['thetao'],
                         gridnow['p'].broadcast_like(dsnow['so']),
                         gridnow['geolon'],
                         gridnow['geolat'],
                         verbose=True)
    
    g = dsnow['so'].copy(data=g_vals).assign_attrs(g_attrs).transpose()
    g.name='gamma_n'
    g.to_netcdf(outpath+'gamma_n.nc')
    
print('---       FIN        ---')