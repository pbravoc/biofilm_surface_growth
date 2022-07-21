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
p1 = @df df boxplot(:strain, :α, outliers=false, group=:strain, ylim=(0.0, 1.5), ylabel="α [1/hr]")
p2 = @df df boxplot(:strain, :β, outliers=false, group=:strain, ylim=(0.0, 0.12), ylabel="β [1/hr]")
p3 = @df df boxplot(:strain, :L, outliers=false, group=:strain, ylim=(0.0, 70.0), ylabel="L [μm]" )
p4 = @df df boxplot(:strain, :h_max, outliers=false, group=:strain, ylim=(0.0, 1100.0), ylabel="h_max [μm]")
plot(p1, p2, p3, p4, legend=false, layout=(4,1), grid=false, xrotation=50, size=(400, 750), left_margin=5mm)
savefig("figs/fig5/parameters.svg")
##
tf = filter(x->x.strain .== "gob33", df)
@df tf density(:L)

##
@df df scatter(:α, :β, group=:strain, xlim=(0, 3.0))

## Correlations between parameters
ols_ab = r2(lm(@formula(α_mean  ~ β_mean), gf_summary))
ols_al = r2(lm(@formula(α_mean  ~ L_mean), gf_summary))
ols_bl = r2(lm(@formula(β_mean  ~ L_mean), gf_summary))
p1 = @df gf_summary scatter(:α_mean, :β_mean, group=:strain,
                            #series_ann = text.(:strain, :top),
                            xerror=(:α_mean - :α_ql, :α_qh - :α_mean),
                            yerror=(:β_mean - :β_ql, :β_qh - :β_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="β[μm/hr]",
                            title="α-β (R²=0.269)",
                            xlim=(0.3, 2.0), ylim=(0.005, 0.09))
p2 = @df gf_summary scatter(:α_mean, :L_mean, group=:strain, 
                            xerror=(:α_mean - :α_ql, :α_qh - :α_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="α[μm/hr]", ylabel="L[μm]",
                            title="α-L (R²=0.313)",
                            xlim=(0.3, 1.5), ylim=(2, 60))
p3 = @df gf_summary scatter(:β_mean, :L_mean, group=:strain, 
                            xerror=(:β_mean - :β_ql, :β_qh - :β_mean),
                            yerror=(:L_mean - :L_ql, :L_qh - :L_mean), 
                            legend=false, markersize=5, grid=false,
                            xlabel="β[μm/hr]", ylabel="L[μm]",
                            xlim=(0.005, 0.09), ylim=(2, 60),
                            title="β-L (R²=0.054)")##

p4 = @df df scatter(:α, :β, group=:strain, alpha=0.25,
                    legend=false, markersize=2, grid=false,
                    xlabel="α[μm/hr]", ylabel="β[μm/hr]",
                    xlim=(0.3, 2.0), ylim=(0.005, 0.09),
                    markerstrokecolor=:auto,
                    title="α-β (R²=0.017)")                          
p5 = @df df scatter(:α, :L, group=:strain, alpha=0.25,
                    legend=false, markersize=2, grid=false,
                    xlabel="α[μm/hr]", ylabel="L[μm]",
                    xlim=(0.3, 1.5), ylim=(2, 60),
                    markerstrokecolor=:auto,
                    title="α-L (R²=0.067)")
p6 = @df df scatter(:β, :L, group=:strain, alpha=0.25,
                    legend=false, markersize=2, grid=false,
                    xlabel="β[μm/hr]", ylabel="L[μm]",
                    xlim=(0.005, 0.09), ylim=(2, 60),
                    markerstrokecolor=:auto,
                    title="β-L (R²=0.004)")
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(700, 450), dpi=600)
#savefig("figs/fig5/parameters_correlation_summary.png")

##
ols_ab = r2(lm(@formula(α  ~ β), df))
ols_al = r2(lm(@formula(α  ~ L), df))
ols_bl = r2(lm(@formula(β  ~ L), df))
##
@df df scatter(:α, :β, :L, group=:strain, alpha=1,
                legend=true, markersize=0.3, grid=false,
                xlabel="α[μm/hr]", ylabel="β[μm/hr]", zlabel="L[μm]",
                xlim=(0.3, 2.0), ylim=(0.005, 0.09), zlim=(2, 60),
                markerstrokecolor=:auto, size=(700, 600)) 
savefig("figs/fig5/parameters_correlation_summary.html")
       