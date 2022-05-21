#=
This is the figure that shows all the parameters, h_max 
for all the strains in the same plot. Maybe there's a way to 
join this figure w number 6 (as insets?)
=#

using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, Plots.Measures
using ColorSchemes

df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
##
p1 = @df df boxplot(:strain, :α, outliers=false, group=:strain, ylim=(0.0, 1.5), title="α")
p2 = @df df boxplot(:strain, :β, outliers=false, group=:strain, ylim=(0.0, 0.1), title="β")
p3 = @df df boxplot(:strain, :L, outliers=false, group=:strain, ylim=(0.0, 70.0), title="L" )
p4 = @df df boxplot(:strain, :h_max, outliers=false, group=:strain, ylim=(0.0, 1100.0), title="h_max")
plot(p1, p2, p3, p4, legend=false, grid=false, xrotation=50)
savefig("figs/fig5/parameters.svg")

##
lb(x) = quantile(x, .025)
ub(x) = quantile(x, .975)
gf = groupby(df, :strain)
sf = combine(gf, :α .=> lb, :α .=> ub,
            :β .=> lb, :β .=> ub,
            :L .=> lb, :L .=> ub)
print(sf)
##
@df sf plot(:α_mean)
@df sf plot!(:α_lb)
@df sf plot!(:α_ub)