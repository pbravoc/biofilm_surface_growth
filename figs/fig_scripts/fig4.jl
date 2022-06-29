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
bgt_2p = dh_2p(df, "bgt127")
jt_2p = dh_2p(df, "jt305")
gob_2p = dh_2p(df, "gob33")
bgt_3p = dh_3p(df, "bgt127")
jt_3p = dh_3p(df, "jt305")
gob_3p = dh_3p(df, "gob33")
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

#
strain = "jt305"
tf = df_48[df_48.strain .== strain, :]
p4 = plot([0, 300], [0, 300], color=:gray, linestyle=:dash)
@df tf plot!(:avg_height, :interface_long, color=my_colors[3], linewidth=3)
plot!(xlim=(0, 270), ylim=(0, 270), legend=false, xticks=[0, 270], yticks=[0,270])
annotate!(125, -35, "Data [μm]", 8)
annotate!(-40, 135, Plots.text("Model [μm]", 8, rotation = 90 ))


#
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
savefig("figs/fig4/fig4_new.svg")

##

##
savefig("figs/fig4/fig4_48h.svg")

##

my_bgt = h_change("bgt127")
my_jt = h_change("jt305")
my_gob = h_change("gob33")
p2 = plot(my_bgt[1][1], my_bgt[1][2], color=my_colors[2], label="Aeromonas", linewidth=3)
plot!(my_bgt[2][1], my_bgt[2][2], color=my_colors[2], linestyle=:dash, label=false, linewidth=3)
plot!(my_jt[1][1], my_jt[1][2], color=my_colors[3], label="E coli", linewidth=3)
plot!(my_jt[2][1], my_jt[2][2], color=my_colors[3], linestyle=:dash, label=false, linewidth=3)
plot!(my_gob[1][1], my_gob[1][2], color=my_colors[4], label="Yeast(aa)", linewidth=3)
plot!(my_gob[2][1], my_gob[2][2], color=my_colors[4], linestyle=:dash, label=false, linewidth=3)
plot!(xlabel="Height [μm]", ylabel="ΔHeight [μm/hr]", grid=false, legend=false)
##
l = @layout [a{0.6w} [b{0.8h} 
                      grid(1,3)]] 
plot(p1, p2, p3, p4, p5, layout=l, size=(750, 350), left_margin=3mm, bottom_margin=3mm, grid=false)

##
plot(xlabel="Time [hr]", ylabel="Height [μm]", ylim=(0, 1060), legend=:topleft, grid=false)

@df df_pred[df_pred.strain .== "bgt127", :] plot!(:time, :tall, color=my_colors[2], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "jt305", :] plot!(:time, :tall, color=my_colors[3], linewidth=3, linestyle=:solid, label=false)
@df df_pred[df_pred.strain .== "gob33", :] plot!(:time, :tall, color=my_colors[4], linewidth=3, linestyle=:solid, label=false)


##RMSE over time
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
