using Base: UV_HANDLE_TYPE_MAX
using DataFrames, CSV, NPZ 
using Plots, StatsPlots 

data = npzread("/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/additional_data/Equilibrium distribution/2021-11-09_bgt127dist/Clean/profiles.npy")
##
plot(data', label=false, ylim=(0, 300))
vline!([5000, 8500])
##
h_center = mean(data[:, 5000:8500], dims=2)
histogram(h_center)
##
prediction = DataFrame(CSV.File("data/sims/bootstrap/boot_bgt127.csv")).h_max
histogram(prediction)
##
histogram([prediction, h_center], bins=120:10:280, normalize=true, label=["Bootstrap" "Experimental"])
