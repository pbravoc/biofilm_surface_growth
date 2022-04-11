#%%
import h5py
import numpy as np
from glob import glob 
import matplotlib.pyplot as plt
import os
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from scipy.optimize import curve_fit

def subplane(img):
    """Substract a plane from the image"""
    X = np.nanmean(img, axis=0)
    valid = ~np.isnan(X)
    popt, pcov = curve_fit(straightline, np.arange(len(X))[valid], X[valid])
    planeX = straightline(np.arange(len(X)), popt[0], popt[1])
    img2 = np.zeros_like(img)
    for i in range(len(planeX)):
        img2[:, i] = img[:,i] - planeX[i]
    
    Y = np.nanmean(img2, axis=1)
    valid = ~np.isnan(Y)
    popt, pcov = curve_fit(straightline, np.arange(len(Y))[valid], Y[valid])
    planeY = straightline(np.arange(len(Y)), popt[0], popt[1])
    img3 = np.zeros_like(img2)
    for i in range(len(planeY)):
        img3[i, :] = img2[i,:] - planeY[i]
    return img3

def datx2py(file_name):
    """Loads a .datx into Python, credit goes to gkaplan.
    https://gist.github.com/g-s-k/ccffb1e84df065a690e554f4b40cfd3a"""
    def _group2dict(obj):
        return {k: _decode_h5(v) for k, v in zip(obj.keys(), obj.values())}
    def _struct2dict(obj):
        names = obj.dtype.names
        return [dict(zip(names, _decode_h5(record))) for record in obj]
    def _decode_h5(obj):
        if isinstance(obj, h5py.Group):
            d = _group2dict(obj)
            if len(obj.attrs):
                d['attrs'] = _decode_h5(obj.attrs)
            return d
        elif isinstance(obj, h5py.AttributeManager):
            return _group2dict(obj)
        elif isinstance(obj, h5py.Dataset):
            d = {'attrs': _decode_h5(obj.attrs)}
            try:
                d['vals'] = obj[()]
            except (OSError, TypeError):
                pass
            return d
        elif isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.number) and obj.shape == (1,):
                return obj[0]
            elif obj.dtype == 'object':
                return _decode_h5([_decode_h5(o) for o in obj])
            elif np.issubdtype(obj.dtype, np.void):
                return _decode_h5(_struct2dict(obj))
            else:
                return obj
        elif isinstance(obj, np.void):
            return _decode_h5([_decode_h5(o) for o in obj])
        elif isinstance(obj, bytes):
            return obj.decode()
        elif isinstance(obj, list) or isinstance(obj, tuple):
            if len(obj) == 1:
                return obj[0]
            else:
                return obj
        else:
            return obj
    with h5py.File(file_name, 'r') as f:
        h5data = _decode_h5(f)
    return h5data

def straightline(x, a, b):
    return a*x + b

def get_data(datx_file):
    """Returns the Surface and Intensity data from a single .datx file"""
    myh5 = datx2py(datx_file)                      # File is the string with the location of the file
    zsurf = myh5['Data']['Surface']           # Get the surfaces
    zdata = list(zsurf.values())[0]           # Good for fixing stuff later  
    zsurf = zdata['vals']                     # Get the data from the surface group
    zsurf[zsurf == zdata['attrs']['No Data']] = np.nan  # Write no data as NaNs for compatibility
    zint = myh5['Data']['Intensity']          # Get the intensity group
    zint = list(zint.values())[0]['vals'].astype(float)  # Get the data from the intensity grou[]
    zint[zint>200000] = np.nan                # This fixes the regions left out from stitching
    return zsurf, zint

def getcleansurf(img):
    trad = np.nanmean(img, axis=0)
    suma = np.nansum(img, axis=0)
    counts = np.sum(~np.isnan(img), axis=0)
    counts[counts<40] = 0
    f = suma/counts
    f[np.isinf(f)] = np.nan
    return f

#%%
files = ["/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_001.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_004.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_007.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_010.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_013.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_016.datx",
"/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_019.datx"]
t = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/times.npy")[1,1:9]

#%%
clean_data = []
for file in files:
    zs, zi = get_data(file)
    t, d = 60, 250
    l,r = 10000, 13000
    cr_2d = subplane(zs[t:d,l:r])
    new_img = np.copy(cr_2d)-np.nanmin(cr_2d)
    clean_data.append(new_img)
t = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/times.npy")[1,1:9]

fig = plt.figure(figsize=(12, 7), constrained_layout=True)
spec = fig.add_gridspec(ncols=1, nrows=7)
for i in range(7):
    ax1 = plt.subplot(spec[i])
    plt.axis('off')
    ax1.imshow(clean_data[i], cmap="viridis", interpolation="none", clim=(0.0, 1000.0))
    ax1.set_aspect('equal')
    txt = plt.text(0.01, 0.75, str(int(t[i]*60))+" minutes", horizontalalignment='left', fontsize=18,
     verticalalignment='center', transform=ax1.transAxes)
    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

norm = plt.Normalize(vmin=0.0, vmax=1000.0)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
cbar_ax = fig.add_axes([1.02, 0.02, 0.05, 0.96])
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(labelsize=14) 
cbar.ax.set_ylabel('Z [$n m$]', rotation=90, fontsize=14)
cbar.set_ticks([0, 200,400, 600, 800, 1000])
cmap = plt.cm.ScalarMappable(norm=norm, cmap='viridis')
cmap.set_array([])
plt.savefig("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/figs/fig1/top_view_nanometer.svg")
# %%
