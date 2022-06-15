##
using DataFrames, CSV
using Statistics, NaNMath, StatsBase
using Plots, StatsPlots
using Optim

## This is to do bootstrap fitting on 48h data
strain_name = "sw519"
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== strain_name && x.time .< 48 &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
interface(h, p) = p[1]*min(h, p[3]) - p[2]*h

h_array = Array(0.0:0.1:200.0)
h_change = [interface(h, [0.6, 20.0, 0.05]) for h in h_array]

@df df scatter(:avg_height, :slope)
plot!(h_array, h_change)
function minimize_function(x)
    h_array = df.avg_height
    h_data = df.slope 
    h_predicted = [interface(h, x) for h in h_array]
    loss = sum(abs.(h_data-h_predicted).^2)
    return loss
end

minimize_function([0.6, 20.0, 0.06])
x0 = [0.6, 20.0, 0.06]
res = optimize(minimize_function, x0)
h_optim = [interface(h, Optim.minimizer(res)) for h in h_array]
@df df scatter(:avg_height, :slope, yerror=:slope_error, color=:black)
plot!(h_array, h_optim, ylim=(0, 10), xlim=(0, 150), linewidth=3)
print(Optim.minimizer(res))
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
##
df = DataFrame(CSV.File("data/timelapses/model_predictionsv2.csv"))
df = filter(x-> x.strain .== "bgt127" && x.time .< 48, df)
@df df scatter(:time, :avg_height, yerror=:std_height, 
               color=:gray, alpha=0.5, label="Data")
@df df plot!(:time, [:interface, :interface_alt], color=[1 2], linewidth=3)
plot!(legend=:topleft)