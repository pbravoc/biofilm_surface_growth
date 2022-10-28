using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Plots.Measures
using LsqFit

function forward_difference(y, x)
    L = length(y)
    fwd_difference = zeros(L)
    fwd_difference[1:L-1] = (y[2:L]-y[1:L-1]) ./ (x[2:L]-x[1:L-1])
    fwd_difference[L] = NaN
    return fwd_difference
end

function smooth_height(y, x)
    model(x, p) = p[1] .+ p[2]*x # Linear model
    y_smooth = zeros(length(y))
    slope_mean = zeros(length(y))
    y_smooth .= NaN              # Start them as NaNs and then fill
    slope_mean .= NaN
    dt = 4.0
    for i=1:length(y)
        idx = (x .> x[i]-dt/2) .&& (x .< x[i]+dt/2)
        x_c, y_c = x[idx], y[idx]
        p_guess = [y_c[1], (y_c[end]-y_c[1])/(x_c[end]-x_c[1])]
        fit = curve_fit(model, x_c, y_c, p_guess)
        y_smooth[i] = model(x[i], fit.param)
        slope_mean[i] = fit.param[2]
    end
    return y_smooth, slope_mean
end

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && 
               x.order .% 2 == 0 &&
               x.replicate in ["A"] &&
               !(isnan(x.avg_height)), Df)

df.fwd_base = forward_difference(df.avg_height, df.time)
@df df scatter(:time, :fwd_base)
resolution = 0.2 # μm 
df.cut_height = resolution .* (df.avg_height .÷ resolution)
df.fwd_base = forward_difference(df.avg_height, df.time)
df.fwd_cut = forward_difference(df.cut_height, df.time)

p1 = @df df scatter(:time, [:avg_height, :cut_height], title="Raw Heights", 
                    xlabel="Time [hr]", ylabel="Height [μm]")
p2 = @df df scatter(:avg_height, :fwd_base, title="Forward Difference",
                    xlabel="Height [μm]", ylabel="ΔHeight [μm/hr]")
p2 = @df df scatter!(:cut_height, :fwd_cut)

df.savg_height, df.slope_base = smooth_height(df.avg_height, df.time)
df.scut_height, df.slope_cut = smooth_height(df.cut_height, df.time)
p3 = @df df scatter(:time, [:savg_height, :scut_height], title="Smoothed Heights", 
                    xlabel="Time [hr]", ylabel="Height [μm]")
p4 = @df df scatter(:savg_height, :slope_base, title="Slope",
                    xlabel="Height [μm]", ylabel="ΔHeight [μm/hr]")
p4 = @df df scatter!(:scut_height, :slope_cut)

plot(p1, p2, p3, p4, grid=false, legend=false, size=(600, 500), markersize=3)
@df df scatter(:cut_height, :fwd_cut)
@df df scatter!(:scut_height, :slope_cut)

