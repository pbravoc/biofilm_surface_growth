using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors

pf = DataFrame(CSV.File("data/timelapses/fit_predictions.csv"))

pf = filter(x->x.strain in ["bgt127", "jt305", "gob33"], pf);
#pf = filter(x->x.strain in ["bgt127", "gob33"], pf);
pf.h_std[7] = 0.5*pf.h_std[6]+0.5*pf.h_std[8]

@df pf scatter(:t, :h, yerror=:h_std, grouping=:strain, alpha=0.25,legend=:topleft, markersize=3, label=["Aeromonas" "Yeast (aa)" "E. coli"])
@df pf plot!(:t, :interface_48, group=:strain, color=[1 2 3], linewidth=3, label=false)
@df pf plot!(:t, :interface_all, group=:strain, color=[1 2 3], linewidth=3, linestyle=:dash, label=false)
#@df pf plot!(:t, :logistic_48, group=:strain, color=[1 2 3], linewidth=3, linestyle=:dash, label=false)
#@df pf plot!(:t, :logistic_all, group=:strain, color=[1 2 3], linewidth=3, linestyle=:dash)
plot!([60.0, 70.0], [20.0, 20.0], color=:black, label="Interface (48)")
plot!([60.0, 70.0], [20.0, 20.0], color=:black, linestyle=:dash, label="Interface (all)")

plot!(xlabel="Time [hr]", ylabel="Height [μm]", size=(500, 400), xlim=(-0.5, 48.0), ylim=(-0.5, 265))
#savefig("figs/animations/longtime11.svg")
##
strain, color = "bgt127", 1
p1 = @df pf[pf.strain .== strain,:] plot(:t, :h-:h, ribbon=:h_std, alpha=0.0, legend=:topleft, fillalpha=0.2,color=:black, label=false)
@df pf[pf.strain .== strain,:] plot!(:t, :interface_48-:h, color=color, linewidth=3, label=false)
@df pf[pf.strain .== strain,:] plot!(:t, :interface_all-:h, color=color, linewidth=3, linestyle=:dash, label=false)
##
@df pf[pf.strain .== strain,:] plot(:t, :h-:h, yerror=:h_std, grouping=:strain, alpha=1.0,legend=:topleft, markersize=3)
##
plot(p1, p2, p3, layout=(3,1), size=(400, 500))
savefig("figs/animations/residuals_interface.svg")
##
strain = "jt305"
pf2 = filter(x->x.strain .== strain, pf);
res_interface = mean(sqrt.((pf2.interface_48 .- pf2.h).^2))
res_logistic = mean(sqrt.((pf2.logistic_48 .- pf2.h).^2))
##
pf2 = DataFrame(CSV.File("data/timelapses/longfit_predictions.csv"))
pf2 = filter(x->x.strain in ["bgt127", "jt305", "gob33"], pf2);
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0 && x.replicate in unique(Df.replicate)[4:end]
                && x.strain in ["bgt127", "gob33", "jt305"], Df)
##
@df Df scatter(:time, :avg_height, yerror=:std_height, group=:strain, color=[1 2 3])

@df pf2 plot!(:t, :interface_48, group=:strain, color=[1 2 3], linewidth=3, label=false)
@df pf2 plot!(:t, :interface_all, group=:strain, color=[1 2 3], linewidth=3, label=false, linestyle=:dash)
plot!([-60.0, -70.0], [20.0, 20.0], color=:black, label="Interface (48h)")
plot!([-60.0, -70.0], [20.0, 20.0], color=:black, linestyle=:dash, label="Interface (all)")

plot!(xlabel="Time [hr]", ylabel="Height [μm]", size=(700, 400), xlim=(-0.5, 350.0), ylim=(-0.5, 800), legend=:topleft, left_margin=2mm)
savefig("figs/animations/longtime10.svg")
