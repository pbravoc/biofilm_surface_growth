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
plot!(xlabel="Time [hr]", ylabel="Height [μm]", ylim=(0, 1060), legend=:topleft, grid=false)
plot!([0, 100], [2000, 2000], color=:black, linewidth=3, linestyle=:dash, label="48h")
plot!([0, 100], [2000, 2000], color=:black, linewidth=3, linestyle=:solid, label="All")

@df df_pred[df_pred.strain .== "bgt127", :] plot!(:time, [:tall, :t48], color=my_colors[2], linewidth=3, linestyle=[:solid :dash], label=false)
@df df_pred[df_pred.strain .== "jt305", :] plot!(:time, [:tall, :t48], color=my_colors[3], linewidth=3, linestyle=[:solid :dash], label=false)
@df df_pred[df_pred.strain .== "gob33", :] plot!(:time, [:tall, :t48], color=my_colors[4], linewidth=3, linestyle=[:solid :dash], label=false)

@df df_long[df_long.strain .== "bgt127",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[2],markersize=4, markerstrokewidth=1.5, label="Aeromonas")
@df df_long[df_long.strain .== "jt305",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[3], markersize=4, markerstrokewidth=1.5,marker=:diamond, label="E coli")
@df df_long[df_long.strain .== "gob33",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[4],markersize=4,markerstrokewidth=1.5, marker=:square,label="Yeast (aa)")
boot_df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
boot_df = filter(x->x.strain in ["bgt127", "jt305", "gob33"], boot_df)
plot!(legend=(0.23, 0.95))
annotate!(20, 1000, "A")

p2 = @df boot_df[boot_df.strain.=="bgt127",:] boxplot(["Aeromonas"], :h_max, outliers=false, color=my_colors[2])
@df boot_df[boot_df.strain.=="jt305",:] boxplot!(["E coli"], :h_max, outliers=false, color=my_colors[3])
@df boot_df[boot_df.strain.=="gob33",:] boxplot!(["Yeast (aa)"], :h_max, outliers=false, color=my_colors[4])
plot!(xrotation=45, legend=false, ylim=(0, 1060))
annotate!(1.0, 1000, "B")

l = @layout [
    a{0.8w} b{0.2w}]
plot(p1, p2, layout=l, size=(500, 350), grid=false) #savefig("figs/fig4/fig4_poster.svg")
savefig("figs/fig4/fig4_left.svg")

## Ribbon parameter fit
df_temp = DataFrame(CSV.File("data/sims/bootstrap/time_jt305_unbounded.csv")) # Experimental
df_temp = filter(x->x.β .> 0, df_temp)
df_temp.h_max = df_temp.α .* df_temp.L ./ df_temp.β
ql(x) = quantile(x, [0.25]) # Quantile low
qh(x) = quantile(x, [0.75]) # Quantile high
best_fits = filter(x->x.id .== "best", df_temp)
gf = groupby(df_temp, :time)
gf_summary = combine(gf, [:α=>median, :α=>ql, :α=>qh,
                          :β=>median, :β=>ql, :β=>qh,
                          :L=>median, :L=>ql, :L=>qh,
                          :h_max=>median, :h_max=>ql, :h_max=>qh,
                          :h_sample=>mean])

# Best fits, PT and uncertainty
#
var_α = (gf_summary.α_qh-gf_summary.α_ql)./ best_fits.α
var_β = (gf_summary.β_qh-gf_summary.β_ql)./ best_fits.β
var_L = (gf_summary.L_qh-gf_summary.L_ql)./ best_fits.L
var_h_max = (gf_summary.h_max_qh-gf_summary.h_max_ql)./ best_fits.h_max


##
p_alpha = @df gf_summary plot(:time, :α_ql, fillrange=:α_qh,
                              color=:black, alpha=0.2, ylabel="α [μm/hr]")
@df best_fits plot!(:time, :α , color=:black, linewidth=3, yticks=[0.33, 0.36])
@df gf_summary plot!(twinx(), :time, 100*:α_ql./ best_fits.α[end],
                     fillrange=100*:α_qh./ best_fits.α[end],
                     color=:black, alpha=0.0, ylabel="α %", 
                     yticks=[96, 104])

p_beta = @df gf_summary plot(:time, :β_ql, fillrange=:β_qh,
                             color=:black, alpha=0.2, ylabel="β [μm/hr]")
@df best_fits plot!(:time, :β, color=:black, linewidth=3, yticks=[0.01, 0.03])
plot!()
@df gf_summary plot!(twinx(), :time, 100*:β_ql./ best_fits.β[end],
                     fillrange=100*:β_qh./ best_fits.β[end],
                     color=:black, alpha=0.0, ylabel="β %",
                     yticks=[100, 350])

p_L = @df gf_summary plot(:time, :L_ql, fillrange=:L_qh,
                          color=:black, alpha=0.2, ylabel="L [μm]")
@df best_fits plot!(:time, :L , color=:black, linewidth=3, yticks=[11.0, 17.0])
@df gf_summary plot!(twinx(), :time, 100*:L_ql./ best_fits.L[end],
                     fillrange=100*:L_qh./ best_fits.L[end],
                    color=:black, alpha=0.0, ylabel="L %",
                    yticks=[90, 150])

p_h = @df gf_summary plot(:time, :h_max_ql, fillrange=:h_max_qh, color=:black, alpha=0.2, ylabel="h_max [μm]")
@df best_fits plot!(:time, :h_max, color=:black, linewidth=3, yticks=[150, 450])
@df gf_summary plot!(twinx(), :time, 100*:h_max_ql./ best_fits.h_max[end],
                     fillrange=100*:h_max_qh./ best_fits.h_max[end],
                     color=:black, alpha=0.0, ylabel="h_max %", 
                     yticks=[40, 120])
plot(p_alpha, p_beta, p_L, p_h, legend=false, layout=(4,1), grid=false, size=(350, 500), right_margin=15mm, xlabel="Time [hr]", left_margin=3mm)
savefig("figs/fig4/fig4_right.svg")

##
l = @layout [
    a{0.8w} b{0.2w}]
plot(layout=l)
##
plot(p1, p2, layout=l, size=(500, 350), grid=false) 
#savefig("figs/fig4/fig4_v2.svg")

#@df gf_summary plot(:time, :β_ql, fillrange=:β_qh)
#savefig("figs/fig4/fig4_fit_vals.svg")

##

p3 = @df best_fits plot(:time, :α ./ best_fits.α[end], color=1, linewidth=2, label="α")
@df best_fits plot!(:time, :β ./ best_fits.β[end], color=2, linewidth=2, label="β")
@df best_fits plot!(:time, :L ./ best_fits.L[end], color=3, linewidth=2, label="L")
plot!(ylabel="Accuracy")
@df best_fits plot!(:time, :h_max ./ best_fits.h_max[end],label="h_max", color=4,linewidth=2, grid=false, xticks=[40, 50, 60, 70, 80])

p4 = plot(best_fits.time, [var_α, var_β, var_L], linewidth=2, label=["α" "β" "L"])
plot!(best_fits.time, var_h_max, color=4, linewidth=2, legend=false)
plot!(xlabel="Time [hr]", ylabel="Uncertainty", ylim=(0.0, 1.25), grid=false)
plot(p3, p4, layout=(2,1), size=(350, 350))
plot(p1, p2, p3, p4, layout=l, size=(700, 350),grid=false, left_margin=3mm, bottom_margin=5mm)

##
#p_alpha = @df gf_summary plot(:time, :α_ql, fillrange=:α_qh, color=:black, alpha=0.2, ylabel="α [μm/hr]")
@df best_fits plot(:time, :α , color=:black, linewidth=3)
@df best_fits plot!(twinx(), :time, :α ./ best_fits.α[end] , color=:black, linewidth=3)
plot!(size=(400, 300), right_margin=10mm)
##
## Accuracy uncertainty scatter
p1 = @df best_fits plot(:α ./ best_fits.α[end],var_α,  color=1, marker=:circle, markersize=3)
p1 = @df best_fits plot!( :β ./ best_fits.β[end], var_β,color=2,  marker=:circle, markersize=3)
p1 = @df best_fits plot!(:L ./ best_fits.L[end], var_L, color=3,  marker=:circle, markersize=3)
p1 = @df best_fits plot!(:h_max ./ best_fits.h_max[end], var_h_max, color=4,  marker=:circle, markersize=3)

plot!(xlabel="Accuracy", ylabel="Uncertainty", ylim=(0.0, 1.27), legend=:bottomright,
      size=(350, 300))




## Boxplots parameter fit
df_temp = DataFrame(CSV.File("data/sims/bootstrap/time_jt305_unbounded.csv")) # Experimental
df_temp = filter(x->x.β .> 0, df_temp)
df_temp.h_max = df_temp.α .* df_temp.L ./ df_temp.β
p_alpha = @df df_temp boxplot(:time, :α,  ylim=(0, 0.4), outliers=false, color=ColorSchemes.okabe_ito[1], yticks=[0.0, 0.4], legend=false)
p_beta = @df df_temp[df_temp.time .<66,:] boxplot(:time, :β,  ylim=(0, 0.07), outliers=false, color=:gray, width=0.8)
p_beta = @df df_temp[df_temp.time .>=66,:] boxplot!(:time, :β, outliers=false, color=ColorSchemes.okabe_ito[1], yticks=[0.0, 0.07], legend=false)
p_L = @df df_temp[df_temp.time .<64,:] boxplot(:time, :L,  ylim=(0, 25.0), outliers=false, color=:gray)
p_L = @df df_temp[df_temp.time .>=64,:] boxplot!(:time, :L, outliers=false, color=ColorSchemes.okabe_ito[1], yticks=[0.0, 25.0], legend=false)

p_h = @df df_temp[df_temp.time .<66,:] boxplot(:time, :h_max,  ylim=(0, 750.0), outliers=false, color=:gray, yticks=[0.0, 750], legend=false)
p_h = @df df_temp[df_temp.time .>= 66,:] boxplot!(:time, :h_max,  ylim=(0, 750.0), outliers=false, color=ColorSchemes.okabe_ito[1])

l = @layout [a{0.45w} b{0.1w} grid(4,1)] 
plot(p1, p2, p_alpha, p_beta, p_L, p_h, layout=l, size=(750, 350), left_margin=3mm, bottom_margin=6mm)
#savefig("figs/fig4/fig4_left.svg")

## RMSE over time

strain = "bgt127"
tf = df_48[df_48.strain .== strain, :]
p3 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[2],linestyle=:dash, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[2], linewidth=2 )
plot!(yticks=[-10, 5], legend=false)

strain = "jt305"
tf = df_48[df_48.strain .== strain, :]
p4 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[3],linestyle=:dash, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[3], linewidth=2 )
plot!(yticks=[-6, 6], legend=false)

strain = "gob33"
tf = df_48[df_48.strain .== strain, :]
p5 = @df tf plot(:time, :avg_height-:avg_height, ribbon=:std_height, color=:gray)
@df tf plot!(:time, :interface-:avg_height, color=my_colors[4], linestyle=:dash, linewidth=2 )
@df tf plot!(:time, :interface_long-:avg_height, color=my_colors[4], linewidth=2 )
plot!(yticks=[-10, 10], legend=false)

p_bot = plot(p3, p4, p5, layout=(1,3), size=(700, 80), bottom_margin=3mm, left_margin=3mm)

## RMSE as bar Plots
rmse_data = zeros(3,3)
for i=1:3
    strain = strain_list[i]
    tf = df_48[df_48.strain .== strain, :]
    rmse(y, x) = sqrt(mean((x-y).^2))
    rmse_data[i,1] = rmse(tf.std_height, tf.avg_height .- tf.avg_height)
    rmse_data[i,2] = rmse(tf.interface, tf.avg_height)
    rmse_data[i,3] = rmse(tf.interface_long, tf.avg_height)
end
groupedbar(rmse_data[:,2:3])
##
l = @layout [ a{0.7h} 
    b]
plot(p_top, p_bot, layout=l, size=(600, 300))
##
l = @layout [ [a b]{0.7h} 
    grid(1,3)
]

plot(p1, p2, p3, p4, p5, layout=l, size=(700, 500), bottom_margin=3mm, left_margin=3mm)
#savefig("figs/fig4/fig4_mix.svg")
##
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
##
l = @layout [
    a{0.8w} b] 
plot(p1, p2, layout=l, size=(700, 300), bottom_margin=5mm, grid=false)

## Predicted vs measured
df_long = DataFrame(CSV.File("data/timelapses/longtime_data.csv")) # Experimental
df_predl = DataFrame(CSV.File("data/sims/bootstrap/boot_trajectories_long_ref.csv")) #Long
strain = "bgt127"
tf_ex = filter(x->x.strain .== strain, df_long)
tf_pr = filter(x->x.strain .== strain, df_predl)
pr_error = std(Array(tf_pr[:,1:1000]), dims=2)
scatter(tf_ex.avg_height, tf_pr.t48, 
        xerror=tf_ex.std_height,
        yerror=pr_error,
        color=my_colors[2],markersize=4, 
        markerstrokewidth=1.5, 
        label="Aeromonas")
strain = "jt305"
tf_ex = filter(x->x.strain .== strain, df_long)
tf_pr = filter(x->x.strain .== strain, df_predl)
pr_error = std(Array(tf_pr[:,1:1000]), dims=2)
scatter!(tf_ex.avg_height, tf_pr.t48, 
         xerror=tf_ex.std_height,
         yerror=pr_error,
         color=my_colors[3],markersize=4, 
         marker=:diamond, markerstrokewidth=1.5, 
         label="E coli")
strain = "gob33"
tf_ex = filter(x->x.strain .== strain, df_long)
tf_pr = filter(x->x.strain .== strain, df_predl)
pr_error = std(Array(tf_pr[:,1:1000]), dims=2)
scatter!(tf_ex.avg_height, tf_pr.t48, 
         xerror=tf_ex.std_height,
         yerror=pr_error,
         color=my_colors[4],markersize=4, 
         marker=:square, markerstrokewidth=1.5, 
         label="Yeast (aa)")
plot!([0, 1000], [0, 1000], color=:gray, linewidth=3, linestyle=:dash,
      label=false, legend=:topleft, size=(380, 300))
#savefig("figs/fig4/PM_t48.svg")