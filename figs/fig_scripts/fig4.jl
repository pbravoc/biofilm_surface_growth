#=
Show the confidence intervals for the predictions for 
bgt127, jt305, and gob33. Compare the 48h residuals 
vs the all fit residuals and the small differences.
Maybe plot densities at different points?

TODO:
-
=#
using DataFrames, CSV
using Statistics, NaNMath, LsqFit
using Plots, StatsPlots, Plots.Measures
using ColorSchemes, Colors

interface_dh(h, p) = p[1]*min(h, p[3])-p[2]*h
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
function h_change(strain)
    h_start = 0
    p = Array(filter(x->x.strain .== strain && 
                     x.fit .== "long", pf)[1, 3:5])
    tf = filter(x->x.strain .== strain, df_pred)[:, [:time, :tall]]
    h_middle = tf.tall[97]
    h_end = tf.tall[end]
    h_solid = Array(h_start:0.1:h_middle)
    h_dash = Array(h_middle:0.1:h_end)
    dh_solid =  [interface_dh(h, p) for h in h_solid]
    dh_dash =  [interface_dh(h, p) for h in h_dash]
    return [h_solid, dh_solid], [h_dash, dh_dash]
end

function dh_3p(df, strain)
    tf = filter(x->x.strain .== strain, df)
    heights = reshape(tf.avg_height, (3,7))
    heights_e = reshape(tf.std_height, (3,7))
    times = reshape(tf.time, (3,7))
    dh = zeros(3,7)
    dh_error = zeros(3,7)
    model(x, p) = p[1] .+ p[2]*x # Linear model
    x, y = df.time, df.avg_height
    for l=2:6
        delta_h = []
        for j=1:3
            for i=1:3, k=1:3 # Left, middle, right
                t0 = times[i, l-1]
                x = [times[i, l-1]-t0, times[j, l]-t0, times[k, l+1]-t0]
                y = [heights[i, l-1], heights[j, l], heights[k, l+1]]
                p_guess = [y[1], (y[3]-y[1])/48]
                fit = curve_fit(model, x, y, [0.05,0.9,0.05], p_guess)
                append!(delta_h, fit.param[2])
            end
        dh[j, l] = mean(delta_h)
        dh_error[j,l] = std(delta_h)
        end
    end
    for j=1:3
        l, delta_h = 1, []       # Start
        for k=1:3
            append!(delta_h, (heights[k, l+1]-heights[j, l])/48)
        end
        dh[j, l] = mean(delta_h)
        dh_error[j,l] = std(delta_h)
        l, delta_h = 7, []       # End
        for i=1:3
            append!(delta_h, (heights[j, l]-heights[i, l-1])/48)
        end
        dh[j, l] = mean(delta_h)
        dh_error[j,l] = std(delta_h)
    end
    return heights, heights_e, dh, dh_error
end

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
## A- Timelapse points

alpha_intensity = 0.1
p1 = plot()
plot!(xlabel="Time [hr]", ylabel="Height [μm]", ylim=(0, 1060), legend=:topleft, grid=false)

@df df_pred[df_pred.strain .== "bgt127", :] plot!(:time, :tall, color=my_colors[2], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "jt305", :] plot!(:time, :tall, color=my_colors[3], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "gob33", :] plot!(:time, :tall, color=my_colors[4], linewidth=3, linestyle=:solid, label=false)

@df df_long[df_long.strain .== "bgt127",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[2],markersize=4, markerstrokewidth=1.5, label="Aeromonas")
@df df_long[df_long.strain .== "jt305",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[3], markersize=4, markerstrokewidth=1.5,marker=:diamond, label="E coli")
@df df_long[df_long.strain .== "gob33",:] scatter!(:time, :avg_height, yerror=:std_height, color=my_colors[4],markersize=4,markerstrokewidth=1.5, marker=:square,label="Yeast (aa)")
boot_df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
boot_df = filter(x->x.strain in ["bgt127", "jt305", "gob33"], boot_df)
plot!(legend=(0.23, 0.95), xlim=(0, 336), xticks=Array(0:48:350))
annotate!(20, 1000, "A")

## B- dh v h 
my_bgt = h_change("bgt127")
my_jt = h_change("jt305")
my_gob = h_change("gob33")
bgt_3p = dh_3p(df_long, "bgt127")
jt_3p = dh_3p(df_long, "jt305")
gob_3p = dh_3p(df_long, "gob33")
p2 = plot()
annotate!(880, 11, "B")

hline!([0.0], color=:black, linewidth=2, style=:dash, alpha=0.5)
plot!(my_bgt[1][1], my_bgt[1][2], color=my_colors[2], label="Aeromonas", linewidth=3)
plot!(my_bgt[2][1], my_bgt[2][2], color=my_colors[2], linestyle=:dash, label=false, linewidth=3)
plot!(my_jt[1][1], my_jt[1][2], color=my_colors[3], label="E coli", linewidth=3)
plot!(my_jt[2][1], my_jt[2][2], color=my_colors[3], linestyle=:dash, label=false, linewidth=3)
plot!(my_gob[1][1], my_gob[1][2], color=my_colors[4], label="Yeast(aa)", linewidth=3)
plot!(my_gob[2][1], my_gob[2][2], color=my_colors[4], linestyle=:dash, label=false, linewidth=3)
plot!(xlabel="Height [μm]", ylabel="ΔHeight [μm/hr]", grid=false, legend=false, xlim=(0, 950))
#scatter!(bgt_2p[1], bgt_2p[3], xerror=bgt_2p[2], yerror=bgt_2p[4],color=my_colors[2])
#scatter!(jt_2p[1], jt_2p[3], xerror=jt_2p[2], yerror=jt_2p[4],color=my_colors[3])
#scatter!(gob_2p[1], gob_2p[3], xerror=gob_2p[2],yerror=jt_2p[4],color=my_colors[4])
scatter!(bgt_3p[1], bgt_3p[3], xerror=bgt_3p[2], yerror=bgt_3p[4],color=my_colors[2])
scatter!(jt_3p[1], jt_3p[3], xerror=jt_3p[2], yerror=jt_3p[4],color=my_colors[3], marker=:diamond)
scatter!(gob_3p[1], gob_3p[3], xerror=gob_3p[2],yerror=gob_3p[4],color=my_colors[4], marker=:square)

## 48h behavior - zoom
strain = "bgt127"
tf = df_48[df_48.strain .== strain, :]
p3 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:avg_height, :interface_long, color=my_colors[2], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Data [μm]", 8)
annotate!(-40, 135, Plots.text("Model [μm]", 8, rotation = 90 ))
annotate!(50, 230, "C")
strain = "jt305"
tf = df_48[df_48.strain .== strain, :]
p4 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:avg_height, :interface_long, color=my_colors[3], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Data [μm]", 8)
annotate!(-40, 135, Plots.text("Model [μm]", 8, rotation = 90 ))
strain = "gob33"
tf = df_48[df_48.strain .== strain, :]
p5 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:avg_height, :interface_long, color=my_colors[4], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Data [μm]", 8)
annotate!(-40, 135, Plots.text("Model [μm]", 8, rotation = 90 ))

plot(p3, p4, p5, layout=grid(1,3), grid=false, size=(400, 0.8*150))

##
l = @layout [a{0.6w} [b{0.75h}
                     grid(1,3)] ]
plot(p1, p2, p3, p4, p5, layout=l, size=(750, 350), left_margin=3mm, bottom_margin=3mm, grid=false)
#savefig("figs/fig4/fig4_new.svg")

## Supplemental Figure 1
pf_boot = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
pf_best = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
pf_boot = filter(x->x.β >0.01, pf_boot)
# 
strain = "bgt127"
tf = filter(x->x.strain .== strain, pf_boot)
tf2 = filter(x->x.strain .== strain && x.fit .== "long", pf_best)
p1 = @df tf density(:α, color=:black, linewidth=2, fillcolor=my_colors[2], fill=true, fillalpha=0.5)
vline!([tf2.x1], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 25))
p2 = @df tf density(:β,color=:black, linewidth=2, fillcolor=my_colors[2], fill=true, fillalpha=0.5)
vline!([tf2.x2], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 190))
p3 = @df tf density(:L, color=:black, linewidth=2, fillcolor=my_colors[2], fill=true, fillalpha=0.5)
vline!([tf2.x3], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 1.2))
#   
strain = "jt305"
tf = filter(x->x.strain .== strain, pf_boot)
tf2 = filter(x->x.strain .== strain && x.fit .== "long", pf_best)
p4 = @df tf density(:α, color=:black, linewidth=2, fillcolor=my_colors[3], fill=true, fillalpha=0.5)
vline!([tf2.x1], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 30))
p5 = @df tf density(:β, color=:black, linewidth=2, fillcolor=my_colors[3], fill=true, fillalpha=0.5)
vline!([tf2.x2], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 48))
p6 = @df tf density(:L, color=:black, linewidth=2, fillcolor=my_colors[3], fill=true, fillalpha=0.5, normalized=true)
vline!([tf2.x3], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 0.22))
#
strain = "gob33"
tf = filter(x->x.strain .== strain, pf_boot)
tf2 = filter(x->x.strain .== strain && x.fit .== "long", pf_best)
p7 = @df tf density(:α, color=:black, linewidth=2, fillcolor=my_colors[4], fill=true, fillalpha=0.5)
vline!([tf2.x1], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 90))
p8 = @df tf density(:β, color=:black, linewidth=2, fillcolor=my_colors[4], fill=true, fillalpha=0.5)
vline!([tf2.x2], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 122))
p9 = @df tf density(:L, color=:black, linewidth=2, fillcolor=my_colors[4], fill=true, fillalpha=0.5)
vline!([tf2.x3], color=:black, linewidth=3, linestyle=:dash,
        ylim=(0, 6.11))
#
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, legend=false, grid=false, yticks=[])
## Bootstrap on h_max
pf_boot = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
pf_best = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
pf_boot = filter(x->x.β >0.01, pf_boot)
df_boot = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
df_boot = filter(x->x.strain in ["bgt127", "jt305", "gob33"] && 
                    x.β > 0, pf_boot)
df_best = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
max_heights = df_best.x1 .* df_best.x3 ./ df_best.x2
mh_48 = max_heights[[1,3,2]]
mh_all = max_heights[[37,39,38]]
gf = groupby(df_boot, :strain)
percentage_border = 0.025
ql(x) = round(quantile(x, [percentage_border])[1], digits=3) # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3) # Quantile high
gf_summary = combine(gf, [:h_max=>ql, :h_max=>qh])
gf_summary.h48 = mh_48
gf_summary.range = gf_summary.h_max_qh - gf_summary.h_max_ql
gf_summary.range_percentage = 100*gf_summary.range ./ gf_summary.h48
print(gf_summary)
@df df_boot violin(:strain, :h_max, color=:black, fillalpha=0.2, label="Bootstrap")
scatter!([0.5, 1.5 , 2.5], mh_48, markersize=5, color=:black, marker=:circle, label="48h")
scatter!([0.5, 1.5 , 2.5], mh_all, markersize=5, color=ColorSchemes.okabe_ito[1], marker=:utriangle, label="All")
plot!(ylim=(0, 1000), size=(300, 400), xrotation=45, grid=false, 
      xticks=([.5, 1.5, 2.5],["Aeromonas", "Yeast (aa)", "E. coli"]),
      ylabel="Maximum height [μm]", left_margin=3mm)
#savefig("figs/fig4/hmax_bootstrap.svg")
## 88h timelapse parameters
df = DataFrame(CSV.File("data/sims/bootstrap/time_jt305_unbounded.csv"))
df.h_max = df.α .* df.L ./ df.β
df_best = filter(x->x.id == "best", df)
df = filter(x->x.α .> 0 && x.β > 0 && x.L > 0, df)
gf = groupby(df, :time)
percentage_border = 0.1
ql(x) = round(quantile(x, [percentage_border])[1], digits=3) # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3) # Quantile high
gf = combine(gf, [:α=>ql, :α=>qh,
                          :β=>ql, :β=>qh,
                          :L=>ql, :L=>qh,
                          :h_max=>ql, :h_max=>qh])
rel = 0.01*df_best.α[end]
p1 = @df gf plot(:time, :α_ql, fillrange=:α_qh, 
                 linewidth=0, fillalpha=0.5, color=:gray, label=false,
                 xlabel="Time [hr]", ylabel="α [μm/hr]")
@df gf plot!(twinx(), :time, :α_ql/rel, fillrange=:α_qh/rel, 
             linewidth=0, fillalpha=0.0, color=:gray, xticks=[],
             label=false, ylabel="%", yguidefontrotation=90)
@df df_best plot!(:time, :α, linewidth=3, color=:black, label=false)
plot!(left_margin=3mm, right_margin=20mm)


rel = 0.01*df_best.β[end]
p2 = @df gf plot(:time, :β_ql, fillrange=:β_qh, 
                 linewidth=0, fillalpha=0.5, color=:gray, label=false,
                 xlabel="Time [hr]", ylabel="β [μm/hr]")
@df gf plot!(twinx(), :time, :β_ql/rel, fillrange=:β_qh/rel, 
             linewidth=0, fillalpha=0.0, color=:gray, xticks=[],
             label=false, ylabel="%", yguidefontrotation=90)
@df df_best plot!(:time, :β, linewidth=3, color=:black, label=false)
plot!(left_margin=3mm, right_margin=20mm)

rel = 0.01*df_best.L[end]
p3 = @df gf plot(:time, :L_ql, fillrange=:L_qh, 
                 linewidth=0, fillalpha=0.5, color=:gray, label=false,
                 xlabel="Time [hr]", ylabel="L [μm]")
@df gf plot!(twinx(), :time, :L_ql/rel, fillrange=:L_qh/rel, 
             linewidth=0, fillalpha=0.0, color=:gray, xticks=[],
             label=false, ylabel="%", yguidefontrotation=90)
@df df_best plot!(:time, :L, linewidth=3, color=:black, label=false)
plot!(left_margin=3mm, right_margin=20mm)

rel = 0.01*df_best.h_max[end]
p4 = @df gf plot(:time, :h_max_ql, fillrange=:h_max_qh, 
                 linewidth=0, fillalpha=0.5, color=:gray, label=false,
                 xlabel="Time [hr]", ylabel="h max [μm]")
@df gf plot!(twinx(), :time, :h_max_ql/rel, fillrange=:h_max_qh/rel, 
             linewidth=0, fillalpha=0.0, color=:gray, xticks=[],
             label=false, ylabel="%", yguidefontrotation=90)
@df df_best plot!(:time, :h_max, linewidth=3, color=:black, label=false)
plot!(left_margin=3mm, right_margin=20mm)
plot(p1, p2, p3, p4, layout=(4,1), size=(500, 600), grid=false)
#savefig("figs/fig4/bootstrap_88.pdf")
