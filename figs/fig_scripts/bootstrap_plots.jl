using Plots, StatsPlots
using StatsBase
using DataFrames, CSV
#plotlyjs()

function remove_outliers(Df, n)
    new_frame = DataFrame()
    for strain in unique(Df.strain)
        df = filter(x-> x.strain .== strain, Df)
        df.L_idx = sortperm(df.L)
        df = filter(x-> x.L_idx .> n && x.L_idx .< 1000-n, df)
        df = select!(df, Not(:L_idx))
        append!(new_frame, df)
    end
    return new_frame 
end
data = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
data = filter(x-> x.strain in ["bgt127", "gob33", "jt305"], data)
#data = remove_outliers(data, 25)
##
data = filter(x-> x.α .>= 0.0 && x.α .<=3.0, data)
data = filter(x-> x.β .>= 0.0 && x.β .<=0.1, data)
data = filter(x-> x.L .>= 0.0 && x.L .<=80.0, data)
##

scatter!(a[:,1], a[:,2], color=Int.(c), markerstrokecolor=:auto, markersize=1.5, alpha=0.3)
##
@df data scatter(:α, :β, :L, group=:strain, 
                 xlim=(0.0, 3.0), ylim=(0.0, 0.1), zlim=(0.0, 80.0),
                 markersize=1, markerstrokecolor=:auto, color=[1 2 3 4 5 6], 
                 alpha=1.0, xlabel="Alpha [um/hr]", ylabel="Beta [um/hr]", zlabel="L [um]",
                 size=(700, 500))
savefig("figs/bootstrapping/scatter.html")        
#@df data_sum scatter!(:α_mean, :β_mean, :L_mean, group=:strain, color=[1 2 3 4 5 6]) 

## BGT127
tf = filter(x-> x.strain .== "bgt127", data)
α = sort(tf.α)[26:974]
β = sort(tf.β)[26:974]
L = sort(tf.L)[26:974]
h_max = sort(tf.h_max)[26:974]
p1 = density(α, fill=true, fillalpha=0.5, color=1, label=false)
plot!(xticks=[0.73, 0.82, params_df[1,3]], yticks=[])
p2 =density(β, fill=true, fillalpha=0.5, color=1, label=false)
plot!(xticks=[0.045, 0.058], yticks=[])
p3 =density(L, fill=true, fillalpha=0.5, color=1, label=false)
plot!(xticks=[13, 17], yticks=[])
p4 = density(h_max, fill=true, fillalpha=0.5, color=1, label=false)
plot!(xticks=[210, 240], yticks=[])
plot(p1, p2, p3, p4, layout=(1,4), size=(500, 150))
## GOB33
tf = filter(x-> x.strain .== "gob33", data)
α = sort(tf.α)[26:974]
β = sort(tf.β)[26:974]
L = sort(tf.L)[26:974]
h_max = sort(tf.h_max)[26:974]
p5 = density(α, fill=true, fillalpha=0.5, color=2, label=false)
plot!(xticks=[0.292, 0.318], yticks=[])
p6 =density(β, fill=true, fillalpha=0.5, color=2, label=false)
plot!(xticks=[0.0, 0.04], yticks=[])
p7 =density(L, fill=true, fillalpha=0.5, color=2, label=false)
plot!(xticks=[28, 42], yticks=[])
p8 = density(h_max, fill=true, fillalpha=0.5, color=2, label=false)
plot!(xticks=[220, 1700], yticks=[])
plot(p5, p6, p7, p8, layout=(1,4), size=(500, 150))
##
tf = filter(x-> x.strain .== "jt305", data)
α = sort(tf.α)[26:974]
β = sort(tf.β)[26:974]
L = sort(tf.L)[26:974]
h_max = sort(tf.h_max)[26:974]
p9 = density(α, fill=true, fillalpha=0.5, color=3, label=false)
plot!(xticks=[0.27, 0.33], yticks=[])
p10 =density(β, fill=true, fillalpha=0.5, color=3, label=false)
plot!(xticks=[0.0, 0.05], yticks=[])
p11 =density(L, fill=true, fillalpha=0.5, color=3, label=false)
plot!(xticks=[11, 25], yticks=[])
p12 = density(h_max, fill=true, fillalpha=0.5, color=3, label=false)
plot!(xticks=[110, 520], yticks=[])
plot(p9, p10, p11, p12, layout=(1,4), size=(500, 150))
##
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, grid=false, layout=(3,4 ), size=(750, 350))
#savefig("figs/bootstrapping/bootstrap_overview.svg")        
##
params_df = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))