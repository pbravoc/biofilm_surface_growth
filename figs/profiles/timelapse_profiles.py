#%%
import numpy as np
import matplotlib.pyplot as plt
from glob import glob 

data = profiles_folder = ()
np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/profiles_A.npy")[0:59,:]
t = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/times.npy")[0,0:59]
