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
               x.time .<= 48 && x.replicate =="A", Df)
pf = DataFrame(CSV.File("data/sims/f3a_heights_bounded.csv"))
#pf = DataFrame(CSV.File("data/sims/f3a_heights_unbounded.csv"))

myc = [ColorSchemes.gray1[1], ColorSchemes.okabe_ito[1],
                 ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[3]] #okabe&ito(2002)
##
p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=myc[1],
               legend=:topleft, label=false)
@df pf plot!(:time, :nutrient_n, color=myc[4], linewidth=2, label="Nutrient")
@df pf plot!(:time, :logistic_n, color=myc[3], linewidth=2, label="Logistic (n)")
@df pf plot!(:time, :interface_n, color=myc[2], linewidth=2, label="Interface (n)")
@df pf plot!(:time, :logistic, color=myc[3], linestyle=:dash, linewidth=2, label="Logistic")
@df pf plot!(:time, :interface, color=myc[2], linestyle=:dash, linewidth=2, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Height [μm]", grid=false, size=(400, 300), dpi=300)
#savefig("figs/fig3/a_unboundedfit.svg")
#annotate!(2, 205, "A")

#
rmse(x, y) = sqrt(mean(x.-y).^2)
my_rmse = [rmse(x, pf.data) for x in [pf.nutrient_n, pf.logistic_n, pf.interface_n, pf.logistic, pf.interface]]
bar!(#["Nutrient (n)", "Logistic (n)", "Interface (n)", "Logistic", "Interface"], 
     my_rmse, xticks=[],
     label=false, color=[myc[4], myc[3], myc[2], myc[3], myc[2]],
     fillstyle=[nothing, nothing, nothing, :x, :x],
     yscale=:log10, xrotation=40, yguidefontsize=10,
     inset = (1, bbox(0.535, 0.72, 0.43, 0.25)),subplot=2, title="RMSE [μm]",
     titlefontsize=10, grid=false)
plot!(size=(400, 300))
#savefig("figs/fig3/fig3_unbounded.pdf")

##
p2 = @df pf plot(:time, :data - :data, color=myc[1], linestyle=:dash, ribbon=:data_error, fillalpha=0.2, label=false)

@df pf plot!(:time, :nutrient_n - :data, color=myc[4], linewidth=2, label="Nutrient_n")
@df pf plot!(:time, :logistic_n- :data, color=myc[3], linewidth=2, label="Logistic_n")
@df pf plot!(:time, :interface_n- :data, color=myc[2], linewidth=2, label="Interface_n")
@df pf plot!(:time, :logistic- :data, color=myc[3], linestyle=:dash, linewidth=2, label="Nutrient")
@df pf plot!(:time, :interface- :data, color=myc[2], linestyle=:dash, linewidth=2, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Residual [μm]", legend=false)

annotate!(45, 25, "B")
##
mean(abs.(pf.interface_n - pf.data))
##
#plot!(inset = (1, bbox(0.35, 0.6, 0.3, 0.35)), subplot=2)
#@df pf scatter!(:data, :nutrient_n, color=3, markersize=2, markerstrokecolor=:auto, label="nutrient_n", subplot=2)
#@df pf scatter!(:data, :logistic, color=2, linestyle=:dash, markersize=2, markerstrokecolor=:auto, label="Logistic", subplot=2)
#@df pf scatter!(:data, :interface, color=1, linestyle=:dash, markersize=2, markerstrokecolor=:auto, label="Interface", subplot=2)
#@df pf plot!([0.0, 200.0], [0.0, 200.0], color=:black, alpha=0.5, linestyle=:dash, linewidth=2,legend=false, subplot=2)
#plot!(xticks=[], yticks=[], subplot=2)
#
function plot_slope(y, x, dt, c, l)
    smooth, slope, slope_error = smooth_heights(y, x, dt)
    plot!(y, slope, ribbon=slope_error, color=c, label=l, linewidth=2)
end

dt = 4.0
h, dh, dh_e = smooth_heights(pf.data, pf.time, dt)
p3 = plot(xlabel="Height [μm]", ylabel="Δ Height [μm/hr]")
scatter!(h, dh, xerror=pf.data_error, yerror=dh_e, color=:black, alpha=0.75,
         markersize=2, label=false)
#vline!([h[findmax(dh)[2]]], color=:black, linestyle=:dash, label=false)
#hline!([0.0], color=:black, linestyle=:dash, label=false, legend=:right)
plot_slope(pf.interface, pf.time, dt, myc[2], "I")
plot_slope(pf.nutrient_n, pf.time, dt, myc[4], "N (n)")
plot_slope(pf.logistic_n, pf.time, dt, myc[3], "L (n)")
#plot!([500.0, 510.0], [[1,1], [2,2], [3,3]], color=[1 2 3], label=["a" "b" "c"])
plot!(xlim=(-1, 220.0), ylim=(-0.1, 13.0), legend=false)
#savefiyg("figs/figs_temp/fig3_c.svg")
annotate!(200, 11.5, "C")
##
l = @layout [
    a{0.6w} [b{0.5h}  
             c{0.5h}] 
]
plot(p1, p2, p3, size=(700, 350), layout=l, bottom_margin=4mm, left_margin=3mm, grid=false)
#savefig("figs/fig3/fig3.pdf")

##


df = DataFrame(CSV.File("data/sims/f3a_heights_bounded.csv"))
aic(rss, n, k) = log(rss/(n-k)+2*k/n)
rss(expected, observed) = sum( (observed .-expected).^2)
#
n_params = [0, 0, 0, 4, 4, 5, 2, 3]
RSS = []
AIC = []
for i in 4:8
    print(i)
    rss_current = rss(df[:, i], df[:,2])
    aic_current = aic(rss_current, 100, n_params[i])
    print(aic_current)
    append!(RSS, rss_current)
    append!(AIC, aic_current)
end
relative_aic = exp.(AIC .-AIC[5])

#
p1 = bar(["Nutrient", "Logistic (n)", "Interface (n)", "Logistic", "Interface"], AIC, color=[1,2,3,2,3],
    fillstyle=[nothing, nothing, nothing, :/, :/], ylabel="AIC", label=false, 
    size=(300, 300), xrotation=40,bottom_margin=4mm)
p2 = bar(["Nutrient", "Logistic (n)", "Interface (n)", "Logistic", "Interface"], relative_aic, color=[1,2,3,2,3],
    fillstyle=[nothing, nothing, nothing, :/, :/], ylabel="Relative likelihood", label=false, 
    size=(300, 300), xrotation=40,bottom_margin=4mm, yscale=:log10)
plot(p1, p2, size=(500, 250), bottom_margin=10mm, dpi=500)
#savefig("figs/fig3/aic.png")
##
##
arss_interface = rss(tf.interface, tf.avg_height)
rss_logistic = rss(tf.logistic, tf.avg_height)
aic_interface = aic(rss_interface, 99, 3)
aic_logistic = aic(rss_logistic, 99, 2)
##
groupedbar(["Aeromonas", "E coli", "Yeast (aa)"], 
           P_rel, bar_position = :dodge, color=myc', 
           label=["Int_long" "Int_48" "Log_long" "Log_48"],
           alpha=[1.0, 0.5, 1.0, 0.5]', bar_width=0.7,
           ylabel="Height % error [μm]", legend=:topleft)

#plot!([0.1, 0.9], [h_bgt127, h_bgt127], color=:black, linestyle=:dash, label=false)
#plot!([1.1, 1.9], [h_jt305, h_jt305], color=:black, linestyle=:dash, label=false)
#plot!([2.1, 2.9], [h_gob33, h_gob33], color=:black, linestyle=:dash, label=false)
plot!(size=(400, 200), dpi=500)
#savefig("figs/fig5/finalheight_temp_rel.png")
## Poster v 
function plot_slope(y, x, dt, c, l)
    smooth, slope, slope_error = smooth_heights(y, x, dt)
    plot!(y, slope, ribbon=slope_error, color=c, label=l, linewidth=2, subplot=2)
end
p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=myc[1],
               legend=:topleft, label=false, grid=false)
@df pf plot!(:time, :logistic_n, color=myc[3], linewidth=3, label="Logistic (n)")
#@df pf plot!(:time, :interface_n, color=myc[2], linewidth=3, label="Interface (n)")
@df pf plot!(:time, :logistic, color=myc[3], linestyle=:dash, linewidth=3, label="Logistic")
#@df pf plot!(:time, :interface, color=myc[2], linestyle=:dash, linewidth=3, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Height [μm]", legend=false)
p3 = plot!(xlabel="h[μm]", ylabel="Δh [μm/hr]", 
           inset = (1, bbox(0.5, 0.48, 0.48, 0.35)),subplot=2,)
scatter!(h, dh, xerror=pf.data_error, yerror=dh_e, color=:black, alpha=0.75,
         markersize=2, label=false, subplot=2)
#vline!([h[findmax(dh)[2]]], color=:black, linestyle=:dash, label=false)
#hline!([0.0], color=:black, linestyle=:dash, label=false, legend=:right)
#plot_slope(pf.interface, pf.time, dt, myc[2], "I")
#plot_slope(pf.nutrient_n, pf.time, dt, myc[4], "N (n)")
plot_slope(pf.logistic_n, pf.time, dt, myc[3], "L (n)")
#plot!([500.0, 510.0], [[1,1], [2,2], [3,3]], color=[1 2 3], label=["a" "b" "c"])
plot!(xlim=(-1, 220.0), ylim=(-0.1, 13.0), legend=false, subplot=2, grid=false)
plot!(size=(500, 300), background=:transparent)
savefig("figs/fig3/fig3_poster.svg")
