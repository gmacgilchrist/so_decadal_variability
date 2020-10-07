# Load and collect ecco data and save in netcdf
import xarray as xr
from scipy.io import loadmat
import numpy as np
import glob
import pandas as pd

# Specify some parameters
start = '1992-01-01'
end = '2017-12-31'
times = pd.date_range(start,end,freq='1M')
gammas = np.arange(26.05,28.6,0.1)
# Load ecco data
rootdir = '/local/projects/so_decadal_variability/ECCOv4r4/binned_vol_budget_gamma/'
filename = 'binned_vol_budget_month_dGamman_0.1_'
files = glob.glob(rootdir+filename)
ecco_vals = {}
# Loop through time
for t in np.linspace(1,312,312):
    file = rootdir+filename+str(int(t))+'.mat'
    print(file)
    # Load .mat file
    ecco_tmp = loadmat(file)
    # Find individual terms in loaded .mat file
    for key in ecco_tmp.keys():
        if key[0:3]=='bin':
            if t==1: # If it's the first timestep, preallocate arrays
                ecco_vals[key] = np.empty(shape=(26,312))
            # Sum up loaded data and place in array
            # Turn it into wmt calculation by multiplying by rho0 and 
            # dividing by \delta gamma
            ecco_vals[key][:,int(t)-1] = 1024.5*np.sum(ecco_tmp[key],0)/0.1
            
# Put each term into a DataArray and all terms into a Datatset
print('Loading to Dataset')
ecco = xr.Dataset()
for key in ecco_tmp.keys():
    if key[0:3]=='bin':
        ecco[key] = xr.DataArray(ecco_vals[key],
                                 dims=['gamma_n','time'],
                                 coords={'gamma_n':gammas,'time':times})
        
# Save to a netcdf file
ecco.to_netcdf('data/binned_vol_budget_month_dGamman_0.1.nc')