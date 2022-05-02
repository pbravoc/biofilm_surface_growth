#%%
import numpy as np
import matplotlib.pyplot as plt
from glob import glob 
import pandas as pd

norm = plt.Normalize(vmin=0, vmax=48.0)          # This is the range for the colorbar from 0 to the last point in the times
sm = plt.cm.ScalarMappable(cmap='inferno_r', norm=norm) # here you choose the colorscheme                             # I don't need all the ticks
cmap = plt.cm.ScalarMappable(norm=norm, cmap='inferno_r')
cmap.set_array([])

data_folder = "/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/"
tl_folders = glob(data_folder+"20*")[1:]
replicates = ["A", "B", "C"]
#%%
fig =  plt.figure(figsize=(10, 14.1))
spec5 = fig.add_gridspec(ncols=3, nrows=8)
for i in range(8):
    folder = tl_folders[i]
    df = pd.read_csv(folder+"/cleaning.csv")
    z = df.zoom[0]
    for j in range(3):
        repli = replicates[j]
        t = np.load(folder+"/times.npy")
        A = np.load(folder+"/profiles_B.npy")
        ax = fig.add_subplot(spec5[i, j])

        t = df[df.replicate == repli].time 
        data = np.load(folder+"/profiles_"+repli+".npy")
        x = np.arange(0, data.shape[1]) * 0.17362 * 1e-3 * 50 / z
        x -= x[int(len(x)/2)]
        my_colors = [cmap.to_rgba(j) for j in t]                               # Here I get the individual colors that are associated to the measurements
        for k in range(1, data.shape[0]):                                             # Loop over profiles and plot with the corresponding color
            plt.plot(x, data[k,:], c=my_colors[k], alpha=0.5);
        plt.xlabel("Distance [$\mu m$]", fontsize=14)
        plt.ylabel("Height [$\mu m$]", fontsize=14)
        plt.ylim(0, 700)
        plt.xlim(-8, 8)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
plt.savefig("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/figs/profiles/timelapse_profiles.svg")
#%%

plt.figure(figsize=(10, 3))
norm = plt.Normalize(vmin=0, vmax=48.0)          # This is the range for the colorbar from 0 to the last point in the times
sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=norm) # here you choose the colorscheme                             # I don't need all the ticks
cmap = plt.cm.ScalarMappable(norm=norm, cmap='viridis_r')
cmap.set_array([])
folder = tl_folders[2]
df = pd.read_csv(folder+"/cleaning.csv")
z = df.zoom[0]
repli = "C"
strain = "sw520"
t = df[df.replicate == repli].time 
data = np.load(folder+"/profiles_"+repli+".npy")
x = np.arange(0, data.shape[1]) * 0.17362 * 1e-3 * 50 / z
x -= x[int(len(x)/2)]
my_colors = [cmap.to_rgba(j) for j in t]                               # Here I get the individual colors that are associated to the measurements
for k in range(1, data.shape[0]):                                             # Loop over profiles and plot with the corresponding color
    plt.plot(x, data[k,:], c=my_colors[k], alpha=1.0);
plt.xlabel("Distance from center[$mm$]", fontsize=14)
plt.ylabel("Height [$\mu m$]", fontsize=14)
plt.ylim(0, 330)
plt.xlim(-6.5, 6.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
cbar = plt.colorbar(sm)                                      # make colorbar
cbar.ax.tick_params(labelsize=14) 
cbar.ax.set_ylabel('Time [hr]', rotation=90, fontsize=14) # Label and rotation
cbar.set_ticks(np.arange(0, 48, 8))    
plt.text(-6, 280, strain+" "+repli, fontsize=18)

plt.tight_layout()                      # I don't need all the ticks
plt.savefig("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/figs/profiles/"+strain+repli+".svg", dpi=300, facecolor='white')

# %%
plt.plot(x, np.transpose(data[0:40,:]), color='gray')
plt.plot(x, np.transpose(data[40:78,:]), color='gray')

# %%
# %%
y_corr = 6*(x-0)**2
plt.plot(x, data[52])
plt.plot(x, data[53]-50)
# %%
