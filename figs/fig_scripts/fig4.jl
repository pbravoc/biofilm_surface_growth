#=
Show the confidence intervals for the predictions for 
bgt127, jt305, and gob33. Compare the 48h residuals 
vs the all fit residuals and the small differences.
Maybe plot densities at different points?

TODO:
-
=#
using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors

df_pred = DataFrame(CSV.File("data/sims/bootstrap/boot_trajectories.csv"))
df_long = DataFrame(CSV.File("data/timelapses/longtime_data.csv"))

strain = "bgt127"
##