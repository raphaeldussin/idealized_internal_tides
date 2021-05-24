import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def find_topo_contour(ds, depth=2000):
    """ find indices for depth contour """
    plt.figure()
    contour = plt.contour(ds['D'], depth)
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

#--------------------- narrow channel --------------------------
geom = xr.open_dataset('../narrow_channel/ocean_geometry.nc')

topoi, topoj = find_topo_contour(geom)

refl = np.nan * xr.ones_like(geom['D'])
refl_angle = np.nan * xr.ones_like(geom['D'])
refl_pref = np.nan * xr.ones_like(geom['D'])
refl_dbl = np.nan * xr.ones_like(geom['D'])
trans = np.nan * xr.ones_like(geom['D'])

refl[topoj, topoi] = 1.
refl_angle[topoj, topoi] = 0.
refl_pref[topoj, topoi] = 0.
refl_dbl[topoj, topoi] = 0.
trans[topoj, topoi] = 0.

ds = xr.Dataset()
ds['LON'] = geom['lonh']
ds['LAT'] = geom['lath']
ds['refl'] = refl
ds['refl_angle'] = refl_angle
ds['refl_pref'] = refl_pref
ds['refl_dbl'] = refl_dbl
ds['trans'] = trans

ds.to_netcdf('IWcoefs_narrow_channel.nc')

