using Base: _require_from_serialized
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
using Plots, StatsPlots, Plots.Measures
using ColorSchemes, Colors

my_colors = [ColorSchemes.okabe_ito[8], ColorSchemes.okabe_ito[5],
             ColorSchemes.okabe_ito[4], ColorSchemes.okabe_ito[6]]
df_pred = DataFrame(CSV.File("data/sims/bootstrap/boot_trajectories.csv")) #Long
df_long = DataFrame(CSV.File("data/timelapses/longtime_data.csv")) # Experimental
df_48 = DataFrame(CSV.File("data/timelapses/model_predictions.csv"))          #48 hour
strain_list = ["bgt127", "jt305", "gob33"]
# Get quantiles for plotting
qf = DataFrame()
for strain in strain_list
    tf = filter(x->x.strain .== strain, df_pred)
    qtf = DataFrame()
    quan = [quantile(tf[i, 1:1000], 
            [0.025, 0.25, 0.75, 0.975]) for i=1:size(tf)[1]]
    quan = reduce(hcat, quan)
    qtf.time = tf.time 
    qtf.q025 = quan[1,:]
    qtf.q250 = quan[2,:]
    qtf.q750 = quan[3,:]
    qtf.q975 = quan[4,:]
    qtf.strain = repeat([strain], size(tf)[1])
    append!(qf, qtf)
end
alpha_intensity = 0.1
p1 = plot()

strain, myc = "bgt127", my_colors[2]
@df qf[qf.strain .== strain, :] plot!(:time, :q025, fillrange=:q975, color=myc, fillalpha=alpha_intensity, alpha=0.0, label=false)

strain, myc = "jt305", my_colors[3]
@df qf[qf.strain .== strain, :] plot!(:time, :q025, fillrange=:q975, color=myc, fillalpha=alpha_intensity, alpha=0.0, label=false)


strain, myc = "gob33", my_colors[4]
@df qf[qf.strain .== strain, :] plot!(:time, :q025, fillrange=:q975, color=myc, fillalpha=alpha_intensity, alpha=0.0, label=false)

@df df_pred[df_pred.strain .== "bgt127", :] plot!(:time, :tall, color=my_colors[2], linewidth=3, label=false)
@df df_pred[df_pred.strain .== "jt305", :] plot!(:time, :tall, color=my_colors[3], linewidth=3, label=false)
@df df_pred[df_pred.strain .== "gob33", :] plot!(:time, :tall, color=my_colors[4], linewidth=3, label=false)

@df df_long[df_long.strain .== "bgt127",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[2],markersize=3, marker=:diamond, label="Aeromonas")
@df df_long[df_long.strain .== "jt305",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[3], markersize=3, label="E coli")
@df df_long[df_long.strain .== "gob33",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[4],markersize=3, marker=:square,label="Yeast (aa)")

plot!(xlabel="Time [hr]", ylabel="Residual [μm]", legend=:topleft, grid=false)

pf_log = DataFrame(CSV.File("data/timelapses/fit_params_logistic.csv"))
pf_int = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
P = zeros(4,3)
sh = filter(x->x.time .> 320, df_long)
h_bgt127 = mean(sh[sh.strain .== "bgt127",:].avg_height)
h_jt305 = mean(sh[sh.strain .== "jt305",:].avg_height)
h_gob33 = mean(sh[sh.strain .== "gob33",:].avg_height)
h_ref = [h_bgt127, h_jt305, h_gob33]
for i=1:3
    strain = strain_list[i]
    tf = filter(x->x.strain .== strain && x.fit .== "long", pf_int)
    P[1,i] = (tf.x1 .* tf.x3 ./tf.x2)[1]
    tf = filter(x->x.strain .== strain  && x.fit .== "48h", pf_int)
    P[2,i] = (tf.x1 .* tf.x3 ./tf.x2)[1]
    tf = filter(x->x.strain .== strain  && x.fit .== "long", pf_log)
    P[3,i] = tf.x2[1]
    tf = filter(x->x.strain .== strain  && x.fit .== "48h", pf_log)
    P[4,i] = tf.x2[1]
end
P = P'
P_rel = abs.(ones(3,4) - P ./ h_ref)
R = zeros(4,3)
rmse(y1, y2) = sqrt(mean( (y1 .- y2).^2 ))
for i=1:3
    strain = strain_list[i]
    tf = df_48[df_48.strain .== strain, :]
    R[1,i] = rmse(tf.avg_height, tf.interface_long)
    R[2,i] = rmse(tf.avg_height, tf.interface)
    R[3,i] = rmse(tf.avg_height, tf.logistic_long)
    R[4,i] = rmse(tf.avg_height, tf.logistic)
end
R = R'
P_dif = abs.(P .- h_ref)
myc = [ColorSchemes.okabe_ito[1], ColorSchemes.okabe_ito[1],
       ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[2]]
p2 = plot()
scatter!(R[1,:], P_dif[1,:], color=myc, marker=:circle, alpha=[1.0, 0.5, 1.0, 0.5], label=false)
scatter!(R[2,:], P_dif[2,:], color=myc, marker=:diamond, alpha=[1.0, 0.5, 1.0, 0.5], label=false)
scatter!(R[3,:], P_dif[3,:], color=myc, marker=:square, alpha=[1.0, 0.5, 1.0, 0.5], label=false)
scatter!([-100], [-100], marker=:circle, color=myc[1], label="Interface")
scatter!([-100], [-100], marker=:circle, color=myc[3], label="Logistic")
plot!(xlim=(1, 1e2), ylim=(2, 1e3))
plot!(xlabel="RMSE", ylabel="Max(h) error", size=(400, 400), dpi=500)
plot!(xscale=:log, yscale=:log)
#

strain = "bgt127"
tf = df_48[df_48.strain .== strain, :]
p3 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[2],linestyle=:dot, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[2], linewidth=2 )
plot!(yticks=[-10, 5], legend=false)

strain = "jt305"
tf = df_48[df_48.strain .== strain, :]
p4 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[3],linestyle=:dot, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[3], linewidth=2 )
plot!(yticks=[-6, 6], legend=false)

strain = "gob33"
tf = df_48[df_48.strain .== strain, :]
p5 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[4], linestyle=:dot, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[4], linewidth=2 )
plot!(yticks=[-10, 10], legend=false)

l = @layout [
    a{0.65w} [b{0.8h}  
             c{0.05h}
             d{0.05h}
             e{0.05h}] 
]
plot(p1, p2, p3, p4, p5, layout=l, size=(700, 350), bottom_margin=3mm, left_margin=3mm)
#savefig("figs/fig4/fig4_mix.svg")
##
df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
tf = df[df.strain .== "bgt127", :]
p6 = @df tf density(sort(tf.α)[26:974], color=my_colors[2], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "jt305", :]
@df tf density!(sort(tf.α)[26:974], color=my_colors[3], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "gob33", :]
@df tf density!(sort(tf.α)[26:974], color=my_colors[4], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "bgt127", :]
p7 = @df tf density(sort(tf.β)[26:974], color=my_colors[2], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "jt305", :]
@df tf density!(sort(tf.β)[26:974], color=my_colors[3], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "gob33", :]
@df tf density!(sort(tf.β)[26:974], color=my_colors[4], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "bgt127", :]
p8 = @df tf density(sort(tf.L)[26:974], color=my_colors[2], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "jt305", :]
@df tf density!(sort(tf.L)[26:974], color=my_colors[3], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "gob33", :]
@df tf density!(sort(tf.L)[26:974], color=my_colors[4], linewidth=3, fill=0,fillalpha=0.1)
plot(p6, p7, p8, legend=false, grid=false, layout=(1,3), size=(500, 100), xticks=false, yticks=false)
##
tf = df[df.strain .== "bgt127", :]
p6 = @df tf density(sort(tf.h_max)[26:974], color=my_colors[2], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "jt305", :]
@df tf density!(sort(tf.h_max)[26:974], color=my_colors[3], linewidth=3, fill=0,fillalpha=0.1)
tf = df[df.strain .== "gob33", :]
@df tf density!(sort(tf.h_max)[26:974], color=my_colors[4], linewidth=3, fill=0,fillalpha=0.1)