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
using ColorSchemes
##
df = DataFrame(CSV.File("data/timelapses/model_predictions.csv"))

## Grid plot
strain_names = unique(df.strain)
P = plot(layout=(3,3), size=(700, 600))
colloquial = ["Aeromonas", "E coli", "Yeast (aa)", 
              "Yeast", "V cholerae (wt)", "V cholerae (eps-)",
              "Klebsiella", "B cereus", "S aureus"]
for i=1:9
    strain = strain_names[i]
    tf = filter(x->x.strain .== strain, df)
    @df tf scatter!(:time, :avg_height, 
                    xerror=:time_error, yerror=:std_height, 
                    color=ColorSchemes.Greys_9[7]
                    , subplot=i, 
                    markerstrokecolor=ColorSchemes.Greys_9[7]
                    , markersize=2)
    annotate!((0.07, 0.9), text(colloquial[i], ColorSchemes.Greys_9[7], :left, 8), subplot=i)
    rmse_interface = round(sqrt(mean((tf.avg_height - tf.interface).^2)), digits=2)
    rmse_logistic = round(sqrt(mean((tf.avg_height - tf.logistic).^2)), digits=2)
    str_interface = string("RMSE=", rmse_interface, "μm")
    str_logistic = string("RMSE=", rmse_logistic, "μm")
    #@df tf plot!(:time, :interface, linewidth=3, color=ColorSchemes.okabe_ito[1], subplot=i)
    #annotate!((0.07, 0.81), text(str_interface, ColorSchemes.Greys_9[7], :left, 6), subplot=i)
    @df tf plot!(:time, :logistic, linewidth=3, color=ColorSchemes.okabe_ito[2], subplot=i)
    annotate!((0.07, 0.81), text(str_logistic, ColorSchemes.Greys_9[7], :left, 6), subplot=i)
    plot!(xlabel="Time [hr]", ylabel="Height [μm]", subplot=i, legend=false)
end 
plot(P, bottom_margin=3mm, left_margin=5mm, grid=false, background_color=:transparent)
#savefig("figs/fig5/fig5_logistic.pdf")

## RMSE bar plot
RMSE_logistic = []
RMSE_interface = []
for i=1:9
    strain = strain_names[i]
    tf = filter(x->x.strain .== strain, df)
    rmse_interface = round(sqrt(mean((tf.avg_height - tf.interface).^2)), digits=2)
    rmse_logistic = round(sqrt(mean((tf.avg_height - tf.logistic).^2)), digits=2)
    append!(RMSE_interface, rmse_interface)
    append!(RMSE_logistic, rmse_logistic)
end 
##
groupedbar(colloquial, [RMSE_interface RMSE_logistic], xrotation=45, 
           color=[ColorSchemes.okabe_ito[1] ColorSchemes.okabe_ito[2]],
           label=["Interface" "Logistic"], ylabel="RMSE [μm]", grid=false,
           size=(400,300), legend=:topleft)
savefig("figs/fig5/fig5_RMSE.pdf")
