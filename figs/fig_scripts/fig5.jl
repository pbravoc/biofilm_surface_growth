#=
Generate a multipanel figure where the experimental heights are compared 
to the interface model predictions.

TODO:
- Change strain numbers for descriptive (jt305->ecoli).
- Use real estate on bottom right: fit parameters? RMSE?
=#
using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, Plots.Measures

df = DataFrame(CSV.File("data/timelapses/model_predictions.csv"))
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
strain_names = unique(df.strain)
P = plot(layout=(3,3), size=(800, 500))
colloquial = ["Aeromonas", "E coli", "Yeast (aa)", 
              "Yeast", "V cholerae (wt)", "V cholerae (eps-)",
              "Klebsiella", "B cereus", "S aureus"]
for i=1:9
    strain = strain_names[i]
    tf = filter(x->x.strain .== strain, df)
    @df tf scatter!(:time, :avg_height, 
                    xerror=:time_error, yerror=:std_height, 
                    color=:gray, subplot=i, 
                    markerstrokecolor=:gray, markersize=2)
    @df tf plot!(:time, :interface, linewidth=2, color=ColorSchemes.okabe_ito[1], subplot=i)
    annotate!((0.07, 0.9), text(colloquial[i], :gray, :left, 8), subplot=i)
    plot!(xlabel="Time [hr]", ylabel="Height [μm]", subplot=i, legend=false)
end 
plot(P, bottom_margin=3mm, left_margin=5mm, grid=false)
savefig("figs/fig5/fig5.svg")

