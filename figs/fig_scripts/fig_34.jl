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
myc = [ColorSchemes.gray1[1], ColorSchemes.okabe_ito[1],
                 ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[3]] #okabe&ito(2002)

p1 = @df pf scatter(:time, :data, yerr=:data_error, markersize=2, color=:black,
               legend=(0.2, 0.93), label=false, markerstrokecolor=:black, z_order=:front)
@df pf plot!(:time, :interface, color=myc[2], linewidth=2, label="Interface")
@df pf plot!(:time, :logistic_n, color=myc[3], linewidth=3, label="Logistic")
plot!(xlabel="Time [hr]", ylabel="Height [μm]", grid=false, size=(400, 300), dpi=300,
      xlim=(0, 48.0))
annotate!(2, 205, "A")

#savefig("figs/fig3/fig3_unbounded.pdf")
h, dh, dh_e = smooth_heights(pf.data, pf.time, 4.0)
plot!(inset = (1, bbox(0.47, 0.53, 0.5, 0.35)),subplot=2,
      grid=false, legend=false, xticks=[0,200])
scatter!(h, dh, xerror=pf.data_error, yerror=dh_e, color=:black,
         markersize=0.5, label=false, 
         subplot=2)

sm, sl, se = smooth_heights(pf.interface, pf.time, 4.0)
plot!(sm, sl, color= myc[2], linewidth=2, subplot=2)

sm, sl, se = smooth_heights(pf.logistic_n, pf.time, 4.0)
plot!(sm, sl, color= myc[3], linewidth=2,
      subplot=2)
annotate!(20, 13.5, text("Δh [μm/hr]", 8), subplot=2)
annotate!(100, -4.0, text("h [μm]", 8), subplot=2)
##
p2 =@df pf plot(:time, :data - :data, color=myc[1], linestyle=:dash, ribbon=:data_error, fillalpha=0.2, label=false)

@df pf plot!(:time, :logistic_n- :data, color=myc[3], linewidth=2, label="Logistic_n")
@df pf plot!(:time, :interface- :data, color=myc[2], linewidth=2, label="Interface")
plot!(xlabel="Time [hr]", ylabel="Res [μm]", legend=false, grid=false,
      yticks=[-20, 10], xlim=(0, 48.0))

annotate!(45, 25, "B")
##
interface_dh(h, p) = p[1]*min(h, p[3])-p[2]*h
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))

my_colors = [ColorSchemes.okabe_ito[8], ColorSchemes.okabe_ito[5],
             ColorSchemes.okabe_ito[4], ColorSchemes.okabe_ito[6]]
df_pred = DataFrame(CSV.File("data/sims/bootstrap/boot_trajectories.csv")) #Long
df_long = DataFrame(CSV.File("data/timelapses/longtime_data.csv")) # Experimental
df_48 = DataFrame(CSV.File("data/timelapses/model_predictions.csv"))          #48 hour
strain_list = ["bgt127", "jt305", "gob33"]

# A- Timelapse points

alpha_intensity = 0.1
p3 = plot(xlabel="Time [hr]", ylabel="Height [μm]", ylim=(0, 1060), legend=:topleft, grid=false)

@df df_pred[df_pred.strain .== "bgt127", :] plot!(:time, :tall, color=my_colors[2], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "jt305", :] plot!(:time, :tall, color=my_colors[3], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "gob33", :] plot!(:time, :tall, color=my_colors[4], linewidth=3, linestyle=:solid, label=false)

@df df_long[df_long.strain .== "bgt127",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[2],markersize=3, label="A. veronii")
@df df_long[df_long.strain .== "jt305",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[3], markersize=3, marker=:diamond, label="E. coli")
@df df_long[df_long.strain .== "gob33",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[4],markersize=3,marker=:square,label="S. cerevisiae")
boot_df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
boot_df = filter(x->x.strain in ["bgt127", "jt305", "gob33"], boot_df)
plot!(legend=(0.15, 0.8), xlim=(0, 336), xticks=Array(0:48:350))
annotate!(320, 100, text("A. veronii", :right, 8))
annotate!(320, 390, text("E. coli", :right, 8))
annotate!(320, 950, text("S. cerivisiae", :right, 8))


#annotate!(20, 1000, "B")
##
using GLM 

function my_r2(strain)
    tf = df_48[df_48.strain .== strain, :]
    return round(r2(lm(@formula(avg_height  ~ interface_long), tf)), digits=4)
end

df_48 = DataFrame(CSV.File("data/timelapses/model_predictions.csv"))          #48 hour
strain_list = ["bgt127", "jt305", "gob33"]

# 48h behavior - zoom
strain = "bgt127"
str_r2 = my_r2(strain)
tf = df_48[df_48.strain .== strain, :]
p4 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:interface_long, :avg_height, color=my_colors[2], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Model Height [μm]", 8)
annotate!(-40, 135, Plots.text("Height [μm]", 8, rotation = 90 ))
annotate!(50, 230, "C")
annotate!((0.65, 0.16), text("R²=$str_r2", 7))
strain = "jt305"
str_r2 = my_r2(strain)
tf = df_48[df_48.strain .== strain, :]
p5 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:interface_long, :avg_height, color=my_colors[3], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Model Height [μm]", 8)
annotate!(-40, 135, Plots.text("Height [μm]", 8, rotation = 90 ))
annotate!((0.65, 0.16), text("R²=$str_r2", 7))
strain = "gob33"
str_r2 = my_r2(strain)
tf = df_48[df_48.strain .== strain, :]
p6 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:interface_long, :avg_height, color=my_colors[4], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Model Height [μm]", 8)
annotate!(-40, 135, Plots.text("Height [μm]", 8, rotation = 90 ))
annotate!((0.65, 0.16), text("R²=$str_r2", 7))
plot(p4, p5, p6, layout=grid(1,3), grid=false, size=(400, 0.8*150))
##
l = @layout [[a{0.8h}
               b    ] [c{0.75h}
                     grid(1,3)] ]
plot(p1, p2, p3, p4, p5, p6, layout=l, size=(600, 350))
savefig("fig34_new.svg")
