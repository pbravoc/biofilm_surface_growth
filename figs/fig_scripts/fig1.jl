using Plots, NPZ

## Aeromonas
folder = "2021-06-25_bgt127"
data = npzread("data/timelapses/"*folder*"/profiles_B.npy")[99,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[1,:]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173
plot(x, data, legend=false, color=:black, linewidth=1.5, size=(1050, 300))

# E coli
folder = "2021-08-27_jt305"
x_off = x[end]+3000
data = npzread("data/timelapses/"*folder*"/profiles_B.npy")[84,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[1,:]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

# Yeast (aa)
folder = "2021-09-03_gob33"
x_off = x[end]+2500
data = npzread("data/timelapses/"*folder*"/profiles_B.npy")[102,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[1,:]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

## Yeast
plot()
folder = "2022-01-28_y55"
x_off = 0
data = npzread("data/timelapses/"*folder*"/profiles_A.npy")[142,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173.*5 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

##
# V cholerae (Wt)
folder = "2022-02-11_bh1514"
x_off = x[end]+3000
data = npzread("data/timelapses/"*folder*"/profiles_A.npy")[82,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

# V cholerae (eps-)
folder = "2022-03-23_ea387"
x_off = x[end]+4000
data = npzread("data/timelapses/"*folder*"/profiles_C.npy")[80,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

# Klebsiella
folder = "2022-03-31_cc117"
x_off = x[end]+3000
data = npzread("data/timelapses/"*folder*"/profiles_C.npy")[88,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

# B cereus
folder = "2022-04-21_sw520"
x_off = x[end]+7000
data = npzread("data/timelapses/"*folder*"/profiles_C.npy")[78,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))

# S aureus
folder = "2022-04-29_sw519"
x_off = x[end]+2500
data = npzread("data/timelapses/"*folder*"/profiles_C.npy")[78,:]
t = npzread("data/timelapses/"*folder*"/times.npy")[:,1]
x = (Array(1:length(data)) .- length(data)/2) .* 0.173 .+ x_off
plot!(x, data, legend=false, color=:black, linewidth=1.5, size=(750, 300))


plot!(xticks=[], yticks=[], axis=false, grid=false, background=:transparent)
plot!(size=(2300, 400), aspect_ratio=:equal)
savefig("figs/profiles/profiles_48_equal.svg")