using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes
using LsqFit

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
               x.time .<= 48 && x.replicate =="A", Df)
pf = DataFrame(CSV.File("data/sims/bgt127.csv"))

##
p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=:black,
               legend=:topleft, label="Data")
@df pf plot!(:time, :nutrient_n, color=3, linewidth=2, label="Nutrient_n")
@df pf plot!(:time, :logistic_n, color=2, linewidth=2, label="Logistic_n")
@df pf plot!(:time, :interface_n, color=1, linewidth=2, label="Interface_n")
@df pf plot!(:time, :logistic, color=2, linestyle=:dash, linewidth=2, label="Nutrient")
@df pf plot!(:time, :interface, color=1, linestyle=:dash, linewidth=2, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Heigt [μm]", grid=false)
#savefig("figs/figs_temp/fig3_a.svg")

##
p2 = hline([0.0], color=:black, alpha=0.6, linestyle=:dash)

@df pf plot!(:time, :nutrient_n - :data, color=3, linewidth=2, label="Nutrient_n")
@df pf plot!(:time, :logistic_n- :data, color=2, linewidth=2, label="Logistic_n")
@df pf plot!(:time, :interface_n- :data, color=1, linewidth=2, label="Interface_n")
@df pf plot!(:time, :logistic- :data, color=2, linestyle=:dash, linewidth=2, label="Nutrient")
@df pf plot!(:time, :interface- :data, color=1, linestyle=:dash, linewidth=2, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Residual [μm]")
#savefig("figs/figs_temp/fig3_b.svg")

##
p3 = hline([0.0], color=:black, alpha=0.6, linestyle=:dash)

@df pf plot!(:time, abs.(:nutrient_n - :data)./:data, color=3, linewidth=2, label="Nutrient_n")
@df pf plot!(:time, abs.(:logistic_n- :data)./:data, color=2, linewidth=2, label="Logistic_n")
@df pf plot!(:time, abs.(:interface_n- :data)./:data, color=1, linewidth=2, label="Interface_n")
@df pf plot!(:time, abs.(:logistic- :data)./:data, color=2, linestyle=:dash, linewidth=2, label="Nutrient")
@df pf plot!(:time, abs.(:interface- :data)./:data, color=1, linestyle=:dash, linewidth=2, label="Interface")
plot!(ylim=(0.0, 1.0))
#savefig("figs/figs_temp/fig3_b2.svg")


##
function plot_slope(y, x, dt, c, l)
    smooth, slope, slope_error = smooth_heights(y, x, dt)
    plot!(y, slope, ribbon=slope_error, color=c, label=l)
    #vline!([smooth[findmax(slope)[2]]], color=c, linestyle=:dash, label=false)
    #hline!([maximum(slope)], color=c, linestyle=:dash, label=false)

end

dt = 4.0
h, dh, dh_e = smooth_heights(pf.data, pf.time, dt)
p4 = plot(xlabel="Height [μm]", ylabel="Δ Height [μm/hr]", grid=false)
scatter!(h, dh, xerror=pf.data_error, yerror=dh_e, color=:black, alpha=0.75,
         markersize=3, label="Data ⟨3⟩")
#vline!([h[findmax(dh)[2]]], color=:black, linestyle=:dash, label=false)
#hline!([maximum(dh)], color=:black, linestyle=:dash, label=false, legend=:right)
plot_slope(pf.interface, pf.time, dt, 1, "Interface")
plot_slope(pf.nutrient_n, pf.time, dt, 2, "Nutrients")
plot_slope(pf.logistic_n, pf.time, dt, 3, "Logistic")
#savefig("figs/figs_temp/fig3_c.svg")

##
plot()
@df pf scatter!(:data, :nutrient_n, color=3, linewidth=2, label="Nutrient_n, R²=0.966")
@df pf scatter!(:data, :logistic, color=2, linestyle=:dash, linewidth=2, label="Logistic, R²=0.9199")
@df pf scatter!(:data, :interface, color=1, linestyle=:dash, linewidth=2, label="Interface, R²=0.999")
plot!([0.0, 200.0], [0.0, 200.0], color=:black, alpha=1.0, linestyle=:dash, legend=:topleft, label=false)
plot!(xlabel="Data [μm]", ylabel="Prediction [μm]")
#savefig("figs/figs_temp/fig3_b3.svg")

##
function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    plot!(t, h_avg, ribbon=h_std, 
             fillalpha=0.6, label=strain_name,
             markersize=3, linewidth=2)
end

df = filter(x->x.replicate in ["A", "B", "C"] && x.time .< 48, Df)
df2 = filter(x->x.replicate in unique(Df.replicate)[4:end], Df)

#@df df scatter(:time, :avg_height, markersize=1, legend=false)
plot(xlabel="Time")
get_average(df, "jt305")
get_average(df, "bgt127")
get_average(df, "gob33")
@df filter(x->x.strain == "jt305", df2) scatter!(:time, :avg_height, yerror=:std_height, color=1, markersize=3, label=false)
@df filter(x->x.strain == "bgt127", df2) scatter!(:time, :avg_height, yerror=:std_height,color=2, markersize=3,  label=false)
@df filter(x->x.strain == "gob33", df2) scatter!(:time, :avg_height, yerror=:std_height,color=3,  markersize=3, label=false)
plot!(legend=:topleft, xlabel="Time [hr]", ylabel="Height [μm]")

jt305_lt = DataFrame(CSV.File("data/sims/jt305_lt.csv"))
@df jt305_lt plot!(:time, [:logistic, :interface], color=1, linestyle=[:dash :solid])

bgt127_lt = DataFrame(CSV.File("data/sims/bgt127_lt.csv"))
@df bgt127_lt plot!(:time, [:logistic, :interface], color=2, linestyle=[:dash :solid])

gob33_lt = DataFrame(CSV.File("data/sims/gob33_lt.csv"))
@df gob33_lt plot!(:time, [:logistic, :interface], color=3, linestyle=[:dash :solid])
plot!(ylim=(0, 600))
#savefig("figs/figs_temp/fig3_d1.svg")


