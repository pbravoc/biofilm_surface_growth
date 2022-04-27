#%%
import h5py
import numpy as np
from glob import glob 
import matplotlib.pyplot as plt
import os
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from scipy.optimize import curve_fit

#%%
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

def cleanpoints(x, y):
    idx = y[~np.isnan(y)]
    return x[idx], y[idx]

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

file = "/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/2021-06-25_bgt127/Raw/bgt127_003.datx"
zs, zi = get_data(file)
t, d = 55, 240
l,r = 8100, 11000
splot = subplane(zs[t:d,l:r])
iplot = zi[t:d, l:r]
cr_2d = subplane(zs[t:d,l:r])
m, b = 3.1, -1100
new_img = np.copy(cr_2d)
x = np.arange(new_img.shape[1])*0.173 -95

for i in range(new_img.shape[1]):
    new_img[:,i] -= m*x[i] + b
cr_1d = np.nanmean(new_img, axis=0)
y = np.nanmean(new_img, axis=0)
data = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/profiles_A.npy")[0:59,:]
t = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/times.npy")[0,0:59]

xs = (np.arange(data.shape[1]) - data.shape[1]/2)*0.173

#%%
my_cmap = 'viridis'
fig = plt.figure(figsize=(12, 7), constrained_layout=True)
spec = fig.add_gridspec(ncols=1, nrows=4, height_ratios=[1, 1, 1, 3])

ax1 = plt.subplot(spec[0])
plt.axis('off')
ax1.imshow(iplot, cmap="gray", interpolation="none")
ax1.set_aspect('equal')
plt.text(0.01, 0.75, 'A', horizontalalignment='center', fontsize=18,
     verticalalignment='center', transform=ax1.transAxes)

ax2 = plt.subplot(spec[1])
plt.axis('off')
ax2.imshow(new_img/1000, cmap=my_cmap, interpolation="none")
ax2.set_aspect('equal')
norm = plt.Normalize(vmin=np.nanmin(new_img/1000), vmax=np.nanmax(new_img/1000))
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
cbar = plt.colorbar(sm, aspect=5, shrink = 0.7)
cbar.ax.tick_params(labelsize=14) 
cbar.ax.set_ylabel('Z [$\mu m$]', rotation=90, fontsize=14)
cbar.set_ticks([0, 2,4])
cmap = plt.cm.ScalarMappable(norm=norm, cmap=my_cmap)
cmap.set_array([])
txt = plt.text(0.01, 0.75, 'B', horizontalalignment='center', fontsize=18,
     verticalalignment='center', transform=ax2.transAxes)
txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

ax3 = plt.subplot(spec[2])
plt.plot(x, cr_1d/1000, color="black")
plt.axhline(0, color='black', linestyle='dashed', alpha=0.1)
plt.xlim(x[0], x[-1])
plt.xlabel("Distance [$\mu m$]", fontsize=14)
plt.ylabel("Height [$\mu m$]", fontsize=14)
plt.text(0.01, 0.85, 'C', horizontalalignment='center', fontsize=18,
     verticalalignment='center', transform=ax3.transAxes)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax4 = plt.subplot(spec[3])
norm = plt.Normalize(vmin=0, vmax=t[-1])
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
cbar = plt.colorbar(sm)
cbar.ax.tick_params(labelsize=14) 
cbar.ax.set_ylabel('Time [hr]', rotation=90, fontsize=14)
cbar.set_ticks([0, 4, 8, 12, 16, 20, 24])
cmap = plt.cm.ScalarMappable(norm=norm, cmap=my_cmap)
cmap.set_array([])
my_colors = [cmap.to_rgba(j) for j in t]
for i in range(1, data.shape[0]):
    ax4.plot(xs, data[i,:], c=my_colors[i]);

plt.xlabel("Distance [$\mu m$]", fontsize=14)
plt.ylabel("Height [$\mu m$]", fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(-2, 165)
plt.text(0.01, 0.95, 'D', horizontalalignment='center', fontsize=18,
     verticalalignment='center', transform=ax4.transAxes)
##
plt.savefig("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/figs/fig1/large_space.svg") # Big distance between panels, final touches on inkscape.
# %%
