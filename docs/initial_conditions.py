#!/usr/bin/env python
# coding: utf-8

# # Initial conditions for idealized test cases

# In[1]:


import xarray as xr
import cmip_basins
import seawater


# Using the temperature and salinity observations from World Ocean Atlas 2013, we're building an average T/S profile:

# In[2]:


url_temp = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/archive/data/0114815/public/temperature/netcdf/decav/1.00/'
url_salt = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/archive/data/0114815/public/salinity/netcdf/decav/1.00/'


# Read the yearly data

# In[3]:


period = '0'  # Annual average
cperiod = str(period).zfill(2)

woa13_t = xr.open_dataset(f'{url_temp}/woa13_decav_t{cperiod}_01.nc',
                          decode_times=False, engine='pydap')
woa13_s = xr.open_dataset(f'{url_salt}/woa13_decav_s{cperiod}_01.nc',
                          decode_times=False, engine='pydap')


# In[4]:


woa13_t


# Compute potential temperature

# In[5]:


p = xr.apply_ufunc(seawater.eos80.pres, woa13_t.depth, woa13_t.t_an,
                   dask='parallelized',
                   output_dtypes=[woa13_t.t_an.dtype])

ptemp = xr.apply_ufunc(seawater.eos80.ptmp, woa13_s.s_an, woa13_t.t_an, p,
                       dask='parallelized',
                       output_dtypes=[woa13_t.t_an.dtype])


# Populate a new dataset with potential temperature and salinity:

# In[6]:


woa13 = xr.Dataset()
woa13['ptemp'] = ptemp
woa13['salt'] = woa13_s['s_an']


# Compute the area of the cells:

# In[7]:


def compute_area_regular_grid(ds, Rearth=6378e3):
    """ compute the cells area on a regular grid """
    import numpy as np

    rfac = 2 * np.pi * Rearth / 360

    up = {"nbounds": 1}
    down = {"nbounds": 0}
    if "time" in ds["lon_bnds"].dims:
        up.update({"time": 0})
        down.update({"time": 0})

    dx1d = rfac * (ds["lon_bnds"].isel(up) - ds["lon_bnds"].isel(down))
    dy1d = rfac * (ds["lat_bnds"].isel(up) - ds["lat_bnds"].isel(down))

    dx2d, dy2d = np.meshgrid(dx1d, dy1d)
    _, lat2d = np.meshgrid(ds["lon"].values, ds["lat"].values)

    dx = dx2d * np.cos(2 * np.pi * lat2d / 360)
    dy = dy2d
    area = dx * dy
    return area


# In[8]:


woa13['area'] = xr.DataArray(compute_area_regular_grid(woa13_t),
                             dims=('lat', 'lon'))


# Create a land sea mask and basin codes:

# In[9]:


woa13['mask'] = xr.where(woa13['ptemp'].isel(depth=0).fillna(-9999.) == -9999., 0, 1).squeeze()


# In[10]:


woa13['mask'].plot()


# In[11]:


woa13['basin'] = cmip_basins.generate_basin_codes(woa13, lon='lon', lat='lat',mask='mask')


# In[12]:


woa13['basin'].plot(cmap='nipy_spectral')


# The resulting dataset is:

# In[13]:


woa13


# ## Average T/S profile for the North Pacific

# Use the basin mask and latitude range to extract North Pacific:

# In[14]:


NPAC = woa13.where(woa13['basin'] == 3).sel(lat=slice(0,65))


# Quick check against full domain:

# In[15]:


p = woa13['area'].where(woa13['basin'] != 0).plot()
p.axes.set_ylim([-90, 90])


# In[16]:


p = NPAC['area'].plot()
p.axes.set_ylim([-90, 90])


# Compute a weighted average over lon and lat. Note that weights cannot contain missing value so we're using array from full domain:

# In[17]:


Tprof = NPAC['ptemp'].weighted(woa13['area']).mean(dim=['lon', 'lat'])


# In[18]:


Sprof = NPAC['salt'].weighted(woa13['area']).mean(dim=['lon', 'lat'])


# In[19]:


Tprof


# Compute densities:

# In[20]:


def dens_wright_eos(T, S, p):
  """
  Equation of state for sea water given by Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
  Units: T[degC],S[PSU],p[Pa]
  Returns density [kg m-3]
  """

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4; b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;   c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  l = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  return (p + p0) / (l + al0*(p+p0))


# In[21]:


p = seawater.eos80.pres(woa13['depth'], Tprof)
dens = dens_wright_eos(Tprof, Sprof, p)


# In[22]:


dens


# build the dataset in layer mode:

# In[23]:


NPAC_profile = xr.Dataset()
NPAC_profile['Layer'] = xr.DataArray(data=dens.squeeze().values[:97], dims=('Layer'),
                                     attrs = dict(long_name = "Layer Target Potential Density",
                                                  units = "kg m-3",
                                                  cartesian_axis = "Z",
                                                  positive = "up"))
NPAC_profile['ptemp'] = xr.DataArray(data=Tprof.squeeze().values[:97], dims=('Layer'),
                                    attrs= dict(long_name = "Potential Temperature",
                                                  units = "deg C",
                                                  cartesian_axis = "Z",
                                                  positive = "up"))
NPAC_profile['salt'] = xr.DataArray(data=Sprof.squeeze().values[:97], dims=('Layer'),
                                    attrs = dict(long_name = "Salinity",
                                                  units = "PSU",
                                                  cartesian_axis = "Z",
                                                  positive = "up"))
NPAC_profile['depth'] = xr.DataArray(data=Tprof['depth'].squeeze().values[:97], dims=('Layer'),
                                     attrs = dict(long_name = "Layer average Depth",
                                                  units = "m",
                                                  cartesian_axis = "Z",
                                                  positive = "up"))


# In[24]:


NPAC_profile


# In[25]:


NPAC_profile.to_netcdf('coordinates_NPAC_profile.nc', format='NETCDF3_CLASSIC')

