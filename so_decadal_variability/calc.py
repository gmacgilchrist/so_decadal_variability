import xarray as xr
# Basic calculations to perform on derived data

def _annual(da):
    annual = da.coarsen(time=12,boundary='trim').mean()
    return annual

def _annualanom(da):
    annual = da.coarsen(time=12,boundary='trim').mean()
    mean = da.mean('time')*xr.ones_like(annual)
    annualanom = annual-mean
    return annualanom

def _anom(da):
    mean = da.mean('time')*xr.ones_like(da)
    anom = da-mean
    return anom
