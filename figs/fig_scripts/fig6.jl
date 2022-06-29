#=
This is the figure that shows all the parameters, h_max 
for all the strains in the same plot. Maybe there's a way to 
join this figure w number 6 (as insets?)
=#

using DataFrames, CSV
using Statistics
using Plots, StatsPlots, Plots.Measures
using ColorSchemes

df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
df = filter(x->x.α .> 0 && x.β .>0 && x.L .>0, df)
gf = groupby(df, :strain)
percentage_border = 0.0
ql(x) = round(quantile(x, [percentage_border])[1], digits=3) # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3) # Quantile high
gf_summary = combine(gf, [:α=>ql, :α=>qh,
                          :β=>ql, :β=>qh,
                          :L=>ql, :L=>qh,
                          :h_max=>ql, :h_max=>qh])
print(gf_summary)
##
p1 = @df df boxplot(:strain, :α, outliers=false, group=:strain, ylim=(0.0, 1.5), title="α [μm/hr]")
p2 = @df df boxplot(:strain, :β, outliers=false, group=:strain, ylim=(0.0, 0.1), title="β [μm/hr]")
p3 = @df df boxplot(:strain, :L, outliers=false, group=:strain, ylim=(0.0, 70.0), title="L [μm]" )
p4 = @df df boxplot(:strain, :h_max, outliers=false, group=:strain, ylim=(0.0, 1100.0), title="h_max [μm]")
plot(p1, p2, p3, p4, legend=false, grid=false, xrotation=50)
#savefig("figs/fig5/parameters.pdf")
