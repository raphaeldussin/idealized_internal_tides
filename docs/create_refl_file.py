import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


def find_topo_contour(ds, depth=2000):
    """ find indices for depth contour """
    plt.figure()
    contour = plt.contour(ds['D'], [depth])
    plt.close()
    for item in contour.allsegs:
        if len(item) !=0:
            # discontinuous contours are stored in multiple arrays
            # so we need to concatenate them
            contourxy = None
            for cont in item:
                if len(cont) != 0:
                    if contourxy is None:
                        contourxy = cont
                    else:
                        contourxy = np.concatenate([contourxy, cont])
    xseg_raw = contourxy[:, 0]
    yseg_raw = contourxy[:, 1]
    nseg = len(xseg_raw)
    iseg_fg = np.floor(xseg_raw).astype('int')
    jseg_fg = np.floor(yseg_raw).astype('int')
    return iseg_fg, jseg_fg


#def add_waveguide(ds, xsource=-190, ysource=20, dx=1., dy=1.):
#    """ add a wave guide for source """
#    isource = np.abs(ds['LON'] - (xsource+0.5*dx)).argmin()
#    jsource = np.abs(ds['LAT'] - (ysource+0.5*dy)).argmin()
#    waveguide = np.zeros(ds['refl'].shape)
#    waveguide += ds['refl'].values
#    waveguide[jsource+1, isource+1] = 1
#    waveguide[jsource+1, isource+0] = 1
#    waveguide[jsource+1, isource-1] = 1
#    waveguide[jsource+0, isource-1] = 1
#    waveguide[jsource-1, isource-1] = 1
#    waveguide[jsource-1, isource+0] = 1
#    waveguide[jsource-1, isource+1] = 1
#    ds['refl'] = xr.DataArray(waveguide, dims=ds['refl'].dims)
#
#    return ds


def create_refl_ridge(ds, geom, rcoef=1., depth=1000):
    """ create refl/transmission coef for mid-ocean ridge """
    topoi, topoj = find_topo_contour(geom, depth=depth)
    refl = np.nan * np.ones(geom['D'].shape)
    refl_angle = -999.9* np.ones(geom['D'].shape)
    refl[topoj, topoi] = rcoef
    for j,i in zip(topoj, topoi):
        if i < 50: # mid-basin
            refl_angle[j, i] = 1 * np.pi / 2
        else:
            refl_angle[j, i] = 3 * np.pi / 2

    refl_pref =0. * np.ones(geom['D'].shape)
    refl_dbl = 0. * np.ones(geom['D'].shape)
    trans = 0. * np.ones(geom['D'].shape)
    
    refl_pref[topoj, topoi] = 1.
    refl_dbl[topoj, topoi] = 0.
    trans[topoj, topoi] = 0.

    ds['refl'] = xr.DataArray(refl, dims=geom['D'].dims)
    ds['refl_angle'] = xr.DataArray(refl_angle, dims=geom['D'].dims)
    ds['refl_pref'] = xr.DataArray(refl_pref, dims=geom['D'].dims)
    ds['refl_dbl'] = xr.DataArray(refl_dbl, dims=geom['D'].dims)
    ds['trans'] = xr.DataArray(trans, dims=geom['D'].dims)
    return ds


def create_refl_walls(ds, geom):
    """ create reflection for E-W walls """

    refl = 0 * np.ones(geom['D'].shape)
    refl_angle = -999.9 * np.ones(geom['D'].shape)
    refl_pref =  0* np.ones(geom['D'].shape)
    refl_dbl = 0 * np.ones(geom['D'].shape)
    trans = 0 * np.ones(geom['D'].shape)

    refl[0,:] = 1.
    refl[-1,:] = 1.
    refl_pref[0,:] = 1.
    refl_pref[-1,:] = 1.
    refl_angle[0,:] = np.pi
    refl_angle[-1,:] = 0.
    refl_dbl[0,:] = 0.
    refl_dbl[-1,:] = 0.
    trans[0,:] = 0.
    trans[-1,:] = 0.

    ds['refl'] = xr.DataArray(refl, dims=geom['D'].dims, attrs = {'_FillValue': 1e+20})
    ds['refl_angle'] = xr.DataArray(refl_angle, dims=geom['D'].dims, attrs = {'_FillValue': 1e+20})
    ds['refl_pref'] = xr.DataArray(refl_pref, dims=geom['D'].dims, attrs = {'_FillValue': 1e+20})
    ds['refl_dbl'] = xr.DataArray(refl_dbl, dims=geom['D'].dims, attrs = {'_FillValue': 1e+20})
    ds['trans'] = xr.DataArray(trans, dims=geom['D'].dims, attrs = {'_FillValue': 1e+20})
    return ds


#--------------------- narrow channel --------------------------------
geom = xr.open_dataset('../narrow_channel/ocean_geometry.nc')

ds = xr.Dataset()
ds['LON'] = geom['lonh']
ds['LAT'] = geom['lath']

ds = create_refl_walls(ds, geom)
ds.to_netcdf('IWcoefs_narrow_channel.nc')


#--------------------- long channel --------------------------------
geom = xr.open_dataset('../long_channel/ocean_geometry.nc')

ds = xr.Dataset()
ds['LON'] = geom['lonh']
ds['LAT'] = geom['lath']

ds = create_refl_walls(ds, geom)
ds.to_netcdf('IWcoefs_long_channel.nc')

#--------------------- narrow channel ridge --------------------------
geom = xr.open_dataset('../narrow_channel_ridge/ocean_geometry.nc')

ds = xr.Dataset()
ds['LON'] = geom['lonh']
ds['LAT'] = geom['lath']

ds = create_refl_ridge(ds, geom, rcoef=1., depth=2000)
#ds = add_waveguide(ds)
ds.to_netcdf('IWcoefs_narrow_channel_ridge.nc')


#--------------------- wide channel -----------------------------------
geom = xr.open_dataset('../wide_channel/ocean_geometry.nc')

ds = xr.Dataset()
ds['LON'] = geom['lonh']
ds['LAT'] = geom['lath']

ds = create_refl_walls(ds, geom)
ds.to_netcdf('IWcoefs_wide_channel.nc')
