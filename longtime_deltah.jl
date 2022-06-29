using DataFrames, CSV
using Statistics, NaNMath, LsqFit
using Plots, StatsPlots, Plots.Measures
using ColorSchemes, Colors

function dh_2p(df, strain)
    tf = filter(x->x.strain .== strain, df)
    heights = reshape(tf.avg_height, (3,7))
    heights_e = reshape(tf.std_height, (3,7))
    times = reshape(tf.time, (3,7))
    dh = zeros(3,7)
    dh_error = zeros(3,7)
    for i=1:3
        for j=1:7
            delta_h = []
            for k=1:3
                if j in Array(1:6)
                    append!(delta_h, (heights[k, j+1]-heights[i, j])/48)
                end
                if j in Array(2:7)
                    append!(delta_h, (heights[i, j]-heights[k, j-1])/48)
                end
            end
            #println(delta_h)
            println(mean(delta_h))
            dh[i, j] = mean(delta_h)
            dh_error[i,j] = std(delta_h)
        end
    end
    return heights, heights_e, dh, dh_error
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

function h_change(strain)
    h_start = 0
    p = Array(filter(x->x.strain .== strain && 
                     x.fit .== "long", pf)[1, 3:5])
    tf = filter(x->x.strain .== strain, df_pred)[:, [:time, :tall]]
    h_middle = tf.tall[80]
    h_end = tf.tall[end]
    h_solid = Array(h_start:0.1:h_middle)
    h_dash = Array(h_middle:0.1:h_end)
    dh_solid =  [interface_dh(h, p) for h in h_solid]
    dh_dash =  [interface_dh(h, p) for h in h_dash]
    return [h_solid, dh_solid], [h_dash, dh_dash]
end

interface_dh(h, p) = p[1]*min(h, p[3])-p[2]*h
##
df = DataFrame(CSV.File("data/timelapses/longtime_data.csv")) # Experimental
df_pred = DataFrame(CSV.File("data/sims/bootstrap/boot_trajectories.csv")) #Long
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))

my_colors = [ColorSchemes.okabe_ito[8], ColorSchemes.okabe_ito[5],
             ColorSchemes.okabe_ito[4], ColorSchemes.okabe_ito[6]]

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

##


jt_3p = dh_3p(df, "jt305")

##
scatter(bgt_2p[1], bgt_2p[3], xerror=bgt_2p[2], yerror=bgt_2p[4], color=1)
scatter!(bgt_3p[1], bgt_3p[b3], xerror=bgt_3p[2], yerror=bgt_3p[4], color=2)
