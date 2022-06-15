using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors
using LsqFit
using Plots.Measures

function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    scatter!(t, h_avg, ribbon=h_std, 
             fillalpha=0.1, alpha=0.8, 
             markersize=3)
end

function smooth_heights(y, x, dt)
    model(x, p) = p[1] .+ p[2]*x # Linear model
    y_smooth = zeros(size(y))
    slope_mean = zeros(size(y))
    slope_error = zeros(size(y))
    y_smooth .= NaN              # Start them as NaNs and then fill
    slope_mean .= NaN
    slope_error .= NaN
    if size(df)[1] > 2           # If we have more points!
        for i=1:length(x)
            idx = (x .> x[i]-dt/2) .&& (x .< x[i]+dt/2)
            x_c, y_c = x[idx], y[idx]
            p_guess = [y_c[1], (y_c[end]-y_c[1])/(x_c[end]-x_c[1])]
            fit = curve_fit(model, x_c, y_c, p_guess)
            y_smooth[i] = model(x[i], fit.param)
            slope_mean[i] = fit.param[2]
            slope_error[i] = sqrt(estimate_covar(fit)[2,2])
        end
    end
    return y_smooth, slope_mean, slope_error
end

# Unbounded parameters fitting and simulation
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && x.replicate in ["A", "B", "C"], Df)
pf = DataFrame(CSV.File("data/sims/f3a_heights_bounded.csv"))
#pf = DataFrame(CSV.File("data/sims/f3a_heights_unbounded.csv"))

myc = [ColorSchemes.gray1[1], ColorSchemes.okabe_ito[1],
                 ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[3]] #okabe&ito(2002)
##
#p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=myc[1],
#               legend=:right, label=false)
@df df scatter(:time, :avg_height)
@df pf plot!(:time, :interface, color=myc[2], linestyle=:dash, linewidth=3, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Height [μm]", grid=false, size=(400, 400), dpi=300)
#savefig("figs/fig3/a_unboundedfit.svg")
