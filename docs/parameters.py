
def betaplane(lat0, Rearth=6.378e6):
    """ write f0 and beta for a beta plane centered on lat0 (deg)"""
    import numpy as np
    omega = (2*np.pi/86400.) # earth rotation rate 2 pi / lenght of day
    deg2rad = 2*np.pi/360
    lat0r = deg2rad * lat0
    f0 = 2 * omega * np.sin(lat0r)
    beta = 2 * omega * np.cos(lat0r) / Rearth
    return f0, beta


f0_narrow, beta_narrow = betaplane(20)
print(f'narrow channel: f0 = {f0_narrow}, beta = {beta_narrow}')

f0_wide, beta_wide = betaplane(40)
print(f'wide channel: f0 = {f0_wide}, beta = {beta_wide}')
