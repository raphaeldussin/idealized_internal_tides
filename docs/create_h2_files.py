import xarray as xr
import numpy as np

# narrow channel
ds = xr.Dataset()
ds['h2'] = xr.DataArray(np.zeros((30,100)), dims=('y', 'x'))
ds.to_netcdf('h2_narrow_chan.nc')

ds = xr.Dataset()
ds['h2'] = xr.DataArray(np.zeros((100,30)), dims=('y', 'x'))
ds.to_netcdf('h2_narrow_chan_R1.nc')
# long channel
ds = xr.Dataset()
ds['h2'] = xr.DataArray(np.zeros((30,800)), dims=('y', 'x'))
ds.to_netcdf('h2_long_chan.nc')

# wide channel
ds = xr.Dataset()
ds['h2'] = xr.DataArray(np.zeros((75,100)), dims=('y', 'x'))
ds.to_netcdf('h2_wide_chan.nc')
