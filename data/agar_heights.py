#%%
import h5py
import numpy as np
from glob import glob 
import matplotlib.pyplot as plt
import os

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
    suma = np.nansum(img, axis=0)
    counts = np.sum(~np.isnan(img), axis=0)
    counts[counts<20] = 0
    f = suma/counts
    f[np.isinf(f)] = np.nan
    return f

root_folder = "/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/"
tl_folders = glob(root_folder+"/20*")
folder = tl_folders[0]+"/"
offnames = [folder+"Clean/offsets_A.npy", folder+"Clean/offsets_B.npy", folder+"Clean/offsets_C.npy"]
Cnames = [folder+"Clean/bounds_A.npy", folder+"Clean/bounds_B.npy", folder+"Clean/bounds_C.npy"]
dnames = [folder+"Clean/displacement_A.npy", folder+"Clean/displacement_B.npy", folder+"Clean/displacement_C.npy"]

files = glob(folder+"Raw/*.datx")                                           # Folder with all the .datx
files.sort()
n_replicates = 3                                                            # Number of timelapse replicates (numbered)
n = int(len(glob(folder+"Raw/*[0-9][0-9][0-9].datx"))/n_replicates)         # Timepoints for each replicate
idx = np.arange(n)*n_replicates    # Timepoints x replicates                # Total timelapse measurements
files_split = [[] for i in range(n_replicates)]                             
for i in idx:
    for j in range(n_replicates):
        files_split[j].append(files[i+j])
#%%
first_pass = []
j = 2                                     # Which of the replicates we'll be working on
for i in range(n):
    s = getcleansurf(get_data(files_split[j][i])[0])# Load profiles as a 1-D array
    first_pass.append(s)
# %%

#%%