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

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && x.replicate =="A", Df)
pf = DataFrame(CSV.File("data/sims/f3a_heights_bounded.csv"))

myc = [ColorSchemes.gray1[1], ColorSchemes.okabe_ito[1],
                 ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[3], 
                 ColorSchemes.okabe_ito[7]] #okabe&ito(2002)
##
p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=myc[1],
               legend=:bottomright, label=false)
@df pf plot!(:time, :interface, color=myc[2], linewidth=2, label="Interface")
@df pf plot!(:time, :logistic_n, color=myc[3], linewidth=2, label="Logistic + Nutrient Depletion")
@df pf plot!(:time, :logistic, color=myc[4], linewidth=2, label="Logistic")
@df pf plot!(:time, :nutrient_n, color=myc[5], linewidth=2, label="Nutrient Depletion")
plot!(xlabel="Time [hr]", ylabel="Height [μm]", grid=false, size=(400, 300), dpi=300)
#savefig("figs/fig3/a_unboundedfit.svg")
#annotate!(2, 205, "A")

##
function plot_slope(y, x, dt, c, )
    smooth, slope, slope_error = smooth_heights(y, x, dt)
    plot!(y, slope, color=c, label=l, linewidth=2)
end

dt = 4.0
h, dh, dh_e = smooth_heights(pf.data, pf.time, dt)
p2 = plot(xlabel="Height [μm]", ylabel="Δ Height [μm/hr]")
scatter!(h, dh, xerror=pf.data_error, yerror=dh_e, color=:black, alpha=0.75,
         markersize=1.5, label=false)
sm, sl, se = smooth_heights(pf.interface, pf.time, dt)
plot!(sm, sl, color= myc[2], linewidth=2)
sm, sl, se = smooth_heights(pf.logistic_n, pf.time, dt)
plot!(sm, sl, color= myc[3], linewidth=2)
sm, sl, se = smooth_heights(pf.logistic, pf.time, dt)
plot!(sm, sl, color= myc[4], linewidth=2)
sm, sl, se = smooth_heights(pf.nutrient_n, pf.time, dt)
plot!(sm, sl, color= myc[5], linewidth=2)
plot!(xlim=(-1, 210.0), ylim=(-0.1, 18.5), legend=false)
#savefiyg("figs/figs_temp/fig3_c.svg")
##
df = DataFrame(CSV.File("data/timelapses/model_predictions_all.csv"))
i_fix = findall(isnan.(df.std_height))
for i in i_fix 
    df.std_height[i] = 0.5*df.std_height[i-1] + 0.5*df.std_height[i+1]
end
strain_names = unique(df.strain)
colloquial = ["A. veronii", "E. coli", "S. cerevisiae (aa)", 
              "S. cerevisiae", "V. cholerae (wt)", "V. cholerae (EPS-)",
              "K. pneumoniae", "B. cereus", "S. aureus"]

RMSE_interface = []
RMSE_logisticnd = []
RMSE_logistic = []
RMSE_nutrient = []
##
for i=1:9
    strain = strain_names[i]
    tf = filter(x->x.strain .== strain, df)
    rmse_interface = round(sqrt(mean((tf.avg_height - tf.interface).^2)), digits=2)
    rmse_logisticnd = round(sqrt(mean((tf.avg_height - tf.logistic_nd).^2)), digits=2)
    rmse_logistic = round(sqrt(mean((tf.avg_height - tf.logistic).^2)), digits=2)
    rmse_nutrient = round(sqrt(mean((tf.avg_height - tf.nutrient).^2)), digits=2)
    append!(RMSE_interface, rmse_interface)
    append!(RMSE_logisticnd, rmse_logisticnd)
    append!(RMSE_logistic, rmse_logistic)
    append!(RMSE_nutrient, rmse_nutrient)
end 
##
p3 = groupedbar(colloquial, [RMSE_interface RMSE_logisticnd RMSE_logistic RMSE_nutrient], xrotation=45, 
           color=myc[2:end]',
        ylabel="RMSE [μm]", grid=false,
           size=(400,300), legend=false, ylim=(0, 45))

##
l = @layout [a{0.3w} b{0.3h}  c{0.4w}] 
plot(p1, p2, p3, size=(900, 600), layout=l, bottom_margin=4mm, left_margin=3mm, grid=false, 
     legendfontsize=5)
savefig("figs/fig3/allmodel_comparison.svg")
