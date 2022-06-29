#%%
import numpy as np
import matplotlib.pyplot as plt

data = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2021-06-25_bgt127/profiles_A.npy")
# %%
plt.figure(figsize=(10, 10))
plt.imshow(data, aspect=200, origin='lower')
# %%
plt.plot(np.transpose(data))
# %%
