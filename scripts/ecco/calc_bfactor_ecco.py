# Calculate the b-factor, for relating neutral and in situ potential
# density, following Iudicone et al. (2008)

import numpy as np
import xarray as xr
from xgcm import Grid
import glob
import matplotlib.pyplot as plt
import gsw

# Load grid data and establish parameters for xgcm grid object
rootdir = '/home/gmacgilchrist/data/ECCO/'
filename = 'ECCO-GRID.nc'
grid = xr.open_dataset(rootdir+filename)
# Define the connectivity between faces
face_connections = {'face':
                    {0: {'X':  ((12, 'Y', False), (3, 'X', False)),
                         'Y':  (None,             (1, 'Y', False))},
                     1: {'X':  ((11, 'Y', False), (4, 'X', False)),
                         'Y':  ((0, 'Y', False),  (2, 'Y', False))},
                     2: {'X':  ((10, 'Y', False), (5, 'X', False)),
                         'Y':  ((1, 'Y', False),  (6, 'X', False))},
                     3: {'X':  ((0, 'X', False),  (9, 'Y', False)),
                         'Y':  (None,             (4, 'Y', False))},
                     4: {'X':  ((1, 'X', False),  (8, 'Y', False)),
                         'Y':  ((3, 'Y', False),  (5, 'Y', False))},
                     5: {'X':  ((2, 'X', False),  (7, 'Y', False)),
                         'Y':  ((4, 'Y', False),  (6, 'Y', False))},
                     6: {'X':  ((2, 'Y', False),  (7, 'X', False)),
                         'Y':  ((5, 'Y', False),  (10, 'X', False))},
                     7: {'X':  ((6, 'X', False),  (8, 'X', False)),
                         'Y':  ((5, 'X', False),  (10, 'Y', False))},
                     8: {'X':  ((7, 'X', False),  (9, 'X', False)),
                         'Y':  ((4, 'X', False),  (11, 'Y', False))},
                     9: {'X':  ((8, 'X', False),  None),
                         'Y':  ((3, 'X', False),  (12, 'Y', False))},
                     10: {'X': ((6, 'Y', False),  (11, 'X', False)),
                          'Y': ((7, 'Y', False),  (2, 'X', False))},
                     11: {'X': ((10, 'X', False), (12, 'X', False)),
                          'Y': ((8, 'Y', False),  (1, 'X', False))},
                     12: {'X': ((11, 'X', False), None),
                          'Y': ((9, 'Y', False),  (0, 'X', False))}}}

# Define vertical metrics as negative, to account for descending coordinate
grid['drW'] = -1 * grid.hFacW * grid.drF #vertical cell size at u point
grid['drS'] = -1 * grid.hFacS * grid.drF #vertical cell size at v point
grid['drC'] = -1 * grid.hFacC * grid.drF #vertical cell size at tracer point

metrics = {
    ('X',): ['dxC', 'dxG'], # X distances
    ('Y',): ['dyC', 'dyG'], # Y distances
    ('Z',): ['drW', 'drS', 'drC'], # Z distances
    ('X', 'Y'): ['rA', 'rAz', 'rAs', 'rAw'] # Areas
}

####
# Define all the functions
# Certaintly I can clean this up substantially
####

def _dhorz(ds,xgrid,var,mask=None):
    '''Take horizontal gradient of scalar field given by var.
    Return dictionary of horizontal gradient in each direction.'''
    
    if mask is not None:
        ds[var] = ds[var].where(ds[mask],np.nan)
    
    gx = xgrid.interp(ds[var], 'X')
    gy = xgrid.interp(ds[var], 'Y', boundary='fill')
    dg = xgrid.diff_2d_vector({'X':gx,'Y':gy},boundary='fill')
    
    D = xgrid.interp_2d_vector({'X':xgrid.get_metric(gx,'X'),'Y':xgrid.get_metric(gy,'Y')},boundary='fill')
    
    return {'X':dg['X']/D['X'],'Y':dg['Y']/D['Y']}

def _dvert(ds,xgrid,var):
    gz = xgrid.interp(ds[var],'Z',boundary='extrapolate')
    dgdz = xgrid.derivative(gz,'Z',boundary='extrapolate')
    return dgdz

def grad(ds,xgrid,var,mask=None):
    dC=_dhorz(ds,xgrid,var,mask)
    dC['Z']=_dvert(ds,xgrid,var)
    return dC

def grad_TSg(ds,xgrid):
    dT = grad(ds,xgrid,var='THETA',mask='maskC')
    dS = grad(ds,xgrid,var='SALT',mask='maskC')
    dg = grad(ds,xgrid,var='GAMMAN',mask='maskC')
    return dT, dS, dg
    
def grad_R(dT,dS,ds):
    alpha = gsw.alpha(ds['SALT'],ds['THETA'],ds['Z'])
    beta = gsw.beta(ds['SALT'],ds['THETA'],ds['Z'])

    rho0 = 1024.5
    dR = {}
    dR['X'] = rho0*(-alpha*dT['X'] + beta*dS['X'])
    dR['Y'] = rho0*(-alpha*dT['Y'] + beta*dS['Y'])
    dR['Z'] = rho0*(-alpha*dT['Z'] + beta*dS['Z'])

def absgrad(dC):
    return xr.ufuncs.sqrt(xr.ufuncs.square(dC['X'])+
                          xr.ufuncs.square(dC['Y'])+
                          xr.ufuncs.square(dC['Z']))

def calc_bfactor(ds,xgrid):
    dT,dS,dg = grad_TSg(ds,xgrid)
    dR = grad_R(dT,dS,ds)
    absgradR = absgrad(dR)
    absgradg = absgrad(dg)
    return absgradg/absgradR


# Writing it as a for loop for ease of saving (consistency with other
# ECCO output)
rootdir = '/data2/project/ECCO4v4/'
years = np.arange(1992,2018)
months = np.arange(1,13)
for year in years:
    for month in months:
        gammain = rootdir+'GAMMAN/'+str(year)+'/GAMMAN_'+str(year)+'_'+str(month).zfill(2)+'.nc'
        thetain = rootdir+'nctiles_monthly/THETA/'+str(year)+'/THETA_'+str(year)+'_'+str(month).zfill(2)+'.nc'
        saltin = rootdir+'nctiles_monthly/SALT/'+str(year)+'/SALT_'+str(year)+'_'+str(month).zfill(2)+'.nc'
        
        gamma = xr.open_dataset(gammain)
        theta = xr.open_dataset(thetain)
        salt = xr.open_dataset(saltin)
        
        ds = xr.merge([gamma,theta,salt,grid]).rename({'tile':'face'})
        # Create the grid object for this dataset
        xgrid = Grid(ds, periodic=False, face_connections=face_connections, metrics=metrics)

        b = calc_bfactor(ds,grid)
        b.name = 'BFACTOR'
        print(b)
        
        bout = rootdir+'BFACTOR/'+str(year)+'/BFACTOR_'+str(year)+'_'+str(month).zfill(2)+'.nc'
        #b.to_netcdf(bout)

