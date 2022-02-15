using Plots
pyplot()
using NPZ
using Plots.Measures
##
data = npzread("data/timelapses/2021-06-25_bgt127/profiles_A.npy")
t = npzread("data/timelapses/2021-06-25_bgt127/times.npy")[1,:]
##
mz = zeros(size(data)[2])
x = (Array(1:size(data)[2]) .- size(data)[2]/2) .* 0.173
plot(xlabel="Space [μm]", ylabel="Height [μm]", xlim=(-2600, 2600))
plot(x, data[1:59,:]', line_z = t[1:59]', linewidth=1.2, 
     label=false, color=:turbo, size=(800, 250), 
     botton_margin=2mm, left_margin=2mm, xlabel="Space [μm]", 
     ylabel="Height [μm]", xlim=(-2600, 2600), grid=false, 
     colorbar_title="Time [hr]", colorbar_ticks=[0,4,8,12,16,20,24])
savefig("/run/media/pablo/T7/bgt127_24h.svg")
###