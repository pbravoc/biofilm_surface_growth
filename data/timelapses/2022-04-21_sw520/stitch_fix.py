#%% On profiles A(t=35) and C(t=29), due to the large lateral size
# and relative roughness of the sample. Zygo Mx (v7.5.1) has some
# issues with stitching, prioritizing local slope. New versions 
# are supposed to address this issue.
# Here, we look at the timepoints *after* the stitching issue
# and correct with a polynomial
# For that we center the homeland, and fit a polynomial
# y = a*x^2 + b*x + c
# This way, we can correct just the homeland

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

cleaning_frame = pd.read_csv("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/cleaning.csv")
# %% A 
data = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/profiles_A_uncorrected.npy")
l,r = 27000, 44500
plt.figure()
plt.subplot(211)
plt.imshow(data, aspect=200)
plt.axvline(l, color='red')
plt.axvline(r, color='red')
x = np.arange(-int((r-l)/2), int((r-l)/2))
best_fit = np.polyfit(x, np.transpose(data[:,l:r]), deg=2)
np.savetxt("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/poly_A_best.csv", np.transpose(best_fit), delimiter=",")
new_fit = np.transpose(np.genfromtxt("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/poly_A_new.csv",delimiter=','))
polys = [np.poly1d(best_fit[:,i])(x) for i in range(data.shape[0])]
polys_new = [np.poly1d(new_fit[:,i])(x) for i in range(data.shape[0])]

plt.figure()
plt.subplot(221)
plt.plot(np.transpose(polys))
plt.subplot(222)
plt.plot(np.transpose(data[:, l:r]))
plt.subplot(223)
plt.plot(np.transpose(polys_new))
fix_homeland = np.copy(data[:, l:r])
for i in range(data.shape[0]):
    fix_homeland[i,:] += polys_new[i] - polys[i]
plt.subplot(224)
plt.plot(np.transpose(fix_homeland))
new_data = np.copy(data)
new_data[:, l:r] = fix_homeland
#np.save("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/profiles_A.npy", new_data)

# %% C

data = np.load("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/profiles_C_uncorrected.npy")
l,r = 28000, 45000
plt.figure()
plt.subplot(211)
plt.imshow(data, aspect=200)
plt.axvline(l, color='red')
plt.axvline(r, color='red')
x = np.arange(-int((r-l)/2), int((r-l)/2))
best_fit = np.polyfit(x, np.transpose(data[:,l:r]), deg=2)
np.savetxt("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/poly_C_best.csv", np.transpose(best_fit), delimiter=",")
new_fit = np.transpose(np.genfromtxt("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/poly_C_new.csv",delimiter=','))
polys = [np.poly1d(best_fit[:,i])(x) for i in range(data.shape[0])]
polys_new = [np.poly1d(new_fit[:,i])(x) for i in range(data.shape[0])]
plt.figure()
plt.subplot(221)
plt.plot(np.transpose(polys))
plt.subplot(222)
plt.plot(np.transpose(data[:, l:r]))
plt.subplot(223)
plt.plot(np.transpose(polys_new))
fix_homeland = np.copy(data[:, l:r])
for i in range(data.shape[0]):
    fix_homeland[i,:] += polys_new[i] - polys[i]
plt.subplot(224)
plt.plot(np.transpose(fix_homeland))
new_data = np.copy(data)
new_data[:, l:r] = fix_homeland
np.save("/home/pablo/Documents/Yggdrasil/Files/biofilm_surface_growth/data/timelapses/2022-04-21_sw520/profiles_C.npy", new_data)

# %%
