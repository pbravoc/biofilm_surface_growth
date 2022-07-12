#=
This is the figure that shows all the parameters, h_max 
for all the strains in the same plot. Maybe there's a way to 
join this figure w number 6 (as insets?)
=#
##
using DataFrames, CSV
using Statistics
using Plots, StatsPlots, Plots.Measures
using ColorSchemes

df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
df = filter(x->x.α .> 0 && x.β .>0 && x.L .>0, df)
gf = groupby(df, :strain)
percentage_border = 0.1
ql(x) = round(quantile(x, [percentage_border])[1], digits=3) # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3) # Quantile high
gf_summary = combine(gf, [:α=>ql, :α=>qh, :α=>mean,
                          :β=>ql, :β=>qh, :β=>mean,
                          :L=>ql, :L=>qh, :L=>mean,
                          :h_max=>ql, :h_max=>qh])
print(gf_summary)
##
p1 = @df df boxplot(:strain, :α, outliers=false, group=:strain, ylim=(0.0, 1.5), title="α [μm/hr]")
p2 = @df df boxplot(:strain, :β, outliers=false, group=:strain, ylim=(0.0, 0.1), title="β [μm/hr]")
p3 = @df df boxplot(:strain, :L, outliers=false, group=:strain, ylim=(0.0, 70.0), title="L [μm]" )
p4 = @df df boxplot(:strain, :h_max, outliers=false, group=:strain, ylim=(0.0, 1100.0), title="h_max [μm]")
plot(p1, p2, p3, p4, legend=false, grid=false, xrotation=50)
#savefig("figs/fig5/parameters.pdf")
##
@df df scatter(:α, :β, group=:strain, xlim=(0, 3.0))
##
ols_ab = r2(lm(@formula(α_mean  ~ β_mean), gf_summary))
ols_al = r2(lm(@formula(α_mean  ~ L_mean), gf_summary))
ols_bl = r2(lm(@formula(β_mean  ~ L_mean), gf_summary))

##
p1 = @df gf_summary scatter(:α_mean, :β_mean, group=:strain,
                            #series_ann = text.(:strain, :top),
                            xerror=(:α_mean - :α_ql, :α_qh - :α_mean),
                            yerror=(:β_mean - :β_ql, :β_qh - :β_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="β[μm/hr]",
                            title="α-β (R²=0.269)")
p2 = @df gf_summary scatter(:α_mean, :L_mean, group=:strain, 
                            xerror=(:α_mean - :α_ql, :α_qh - :α_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="L[μm]",
                            title="α-L (R²=0.313)")
p3 = @df gf_summary scatter(:β_mean, :L_mean, group=:strain, 
                            xerror=(:β_mean - :β_ql, :β_qh - :β_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="β[μm/hr]", ylabel="L[μm]",
                            title="β-L (R²=0.054)")##
plot(p1, p2, p3, layout=(1,3), size=(700, 250), dpi=500)
savefig("figs/fig5/parameters_correlation_summary.svg")

##
p1 = @df df scatter(:α, :β, group=:strain,
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="β[μm/hr]")
##                            
p2 = @df gf_summary scatter(:α_mean, :L_mean, group=:strain, 
                            xerror=(:α_mean - :α_ql, :α_qh - :α_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="L[μm]")
p3 = @df gf_summary scatter(:β_mean, :L_mean, group=:strain, 
                            xerror=(:β_mean - :β_ql, :β_qh - :β_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="β[μm/hr]", ylabel="L[μm]")##
plot(p1, p2, p3, layout=(1,3), size=(700, 250))
##
ols_ab = r2(lm(@formula(α  ~ β), df))
ols_al = r2(lm(@formula(α  ~ L), df))
ols_bl = r2(lm(@formula(β  ~ L), df))