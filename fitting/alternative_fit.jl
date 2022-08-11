##
using DataFrames, CSV
using Statistics, NaNMath, StatsBase
using Plots, StatsPlots
using Optim

## This is to do bootstrap fitting on 48h data
strain_name = "bgt127"
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== strain_name && x.time .< 48 &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
interface(h, p) = p[1]*min(h, p[3]) - p[2]*h

h_array = Array(0.0:0.1:200.0)
h_change = [interface(h, [0.6, 0.05, 20.0]) for h in h_array]

@df df scatter(:avg_height, :slope, color=:gray, alpha=0.75, label="Data")
plot!(h_array, h_change, color=ColorSchemes.okabe_ito[1], linewidth=3, label="Interface")
##
function minimize_function(x)
    h_array = df.avg_height
    h_data = df.slope 
    h_predicted = [interface(h, x) for h in h_array]
    loss = sum(abs.(h_data-h_predicted).^2)
    return loss
end

minimize_function([0.6, 20.0, 0.06])
x0 = [0.6, 20.0, 0.06]
res = optimize(minimize_function, x0, store_trace=true, extended_trace=true)
h_optim = [interface(h, Optim.minimizer(res)) for h in h_array]
@df df scatter(:avg_height, :slope, yerror=:slope_error, color=:black)
plot!(h_array, h_optim, ylim=(0, 10), xlim=(0, 150), linewidth=3)
print(Optim.minimizer(res))

##
P = [p.metadata["centroid"] for p in res.trace]
v = [p.value for p in res.trace]
anim = @animate for i=250:1:455
    @df df scatter(:avg_height, :slope, color=:gray, alpha=0.75, label="Data")
    h_optim = [interface(h, P[i]) for h in h_array]
    plot!(h_array, h_optim, color=ColorSchemes.okabe_ito[1], linewidth=3, label="Interface")
    plot!(title=v[i], xlim=(0, 220), ylim=(0, 15), grid=false, xlabel="Height [μm]",
          ylabel="Δ Height [μm/hr]", size=(450, 400), dpi=500)
end
gif(anim, "figs/alt_fitting.gif")

##
strain_list = unique(Df.strain)
P = []
Strain = []
Fit = []
x0 = [0.6, 10.0, 0.06]
for strain in strain_list 
    print(strain)
    df = filter(x-> x.strain .== strain && x.time .<48, Df)              
    res = optimize(minimize_function, x0)
    append!(P, [Optim.minimizer(res)])
    println(Optim.minimizer(res))
    append!(Strain, [strain])
    append!(Fit, ["alt"])
end
##
pf = hcat(DataFrame("strain"=>Strain, "fit"=>Fit),
          DataFrame(Matrix(reduce(hcat, P)'), :auto))
##
CSV.write("data/timelapses/fit_params_alt.csv", pf)
##
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
pf = filter(x-> x.fit in ["48h", "alt"] , pf)
pf.h_max = pf.x1 .* pf.x3 ./ pf.x2
p1 = @df pf groupedbar(:strain, :x1, group=:fit, ylabel="α [μm/hr]")
p2 = @df pf groupedbar(:strain, :x2, group=:fit, legend=false, ylabel="β [μm/hr]")
p3 = @df pf groupedbar(:strain, :x3, group=:fit, legend=false, ylabel="L [μm]")
p4 = @df pf groupedbar(:strain, :h_max, group=:fit, legend=false, ylabel="h_max [μm]")
plot(p1, p2, p3, p4, size=(600, 400), xrotation=45, grid=false)
savefig("figs/alt_fitting_vals.svg")
##
df = DataFrame(CSV.File("data/timelapses/model_predictionsv2.csv"))
df = filter(x-> x.strain .== "sw519" && x.time .< 48, df)
@df df scatter(:time, :avg_height, yerror=:std_height, 
               color=:gray, alpha=0.75, label="Data")
@df df plot!(:time, [:interface, :interface_alt], color=[1 2], linewidth=3, label=["48h" "Alt"])
plot!(legend=:topleft, grid=false, size=(600, 400), xlabel="Time [hr]", ylabel="Height [μm]")
savefig("figs/alt_sw519.svg")
