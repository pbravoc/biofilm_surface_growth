using DataFrames, CSV, NPZ
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Plots.Measures
using SpecialFunctions

G(z, zstar) = z < zstar ? z : zstar
function c(x, t)
    L = sqrt(D*t)
    return erfc.(x / L)
end
function int_c(h)
    L = sqrt(D*dt)
    t1 = L*(-exp.(-h^2 / (L^2)))
    t2 = sqrt(π) * h * erfc(h/L)
    return (t1 + t2 + L)
end

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && x.replicate in ["A", "B", "C"], Df)
x = Array(0.0:0.05:3.0)
D = 0.3
dt = 1.0
y = exp.(-D*x)
myc = c(x, dt)
a = int_c(0.1)
sol_dis = cumsum(c(x, dt)) .- 1
sol_ana = int_c.(x)
sol_dis = G.(x, sqrt(D*dt))
#
my_colors = [colorant"#B3B3B3", colorant"#808080", colorant"#4D4D4D"]'
#
p1 = @df df plot(:time, :avg_height, group=:replicate,
                fillalpha=0.3, ylabel="Height [μm]", xlabel="Time [hr]",
                alpha=1.0,  grid=false, 
                legend=:topleft, color=my_colors,
                marker=:circle, markersize=2.5, markerstrokecolor=:auto,
                linewidth=1)
p1 = hline!([27.5], color=ColorSchemes.okabe_ito[1], linestyle=:dash, linewidth=3, legend=false)
p1 = annotate!(10, 17, text("I", 8, "courier", :left))
p1 = annotate!(10, 38, text("II", 8, "courier", :left))
annotate!(2, 205, "A")
@df df plot!(:time, :smooth_height, group=:replicate,
             alpha=0.8, legend=false, 
             grid=false, color=my_colors,
             marker=:circle, markersize=1, markerstrokecolor=:auto,
             inset = (1, bbox(0.5, 0.45, 0.5, 0.45)), subplot=2,
             linewidth=1.5,yscale=:log10)
#

p2 = @df df plot(:avg_height, :slope, group=:replicate,
            fillalpha=0.3, alpha=0.8,
            ribbon=:slope_error, grid=false, color=my_colors,
            marker=:circle, markersize=3, markerstrokecolor=:auto,
            xlabel="Height [μm]", ylabel="Δ Height [μm/hr]",
            linewidth=2)
p2 = vline!([27.5], color=ColorSchemes.okabe_ito[1], linewidth=3,linestyle=:dash, legend=false)
p2 = annotate!(20, 1, text("I", 8, "courier"))
p2 = annotate!(40, 1, text("II", 8, "courier"))
annotate!(5, 12, "B")
##
@df df plot(:time, :slope, group=:replicate,
            fillalpha=0.3, alpha=0.8,
            ribbon=:slope_error, grid=false,
            marker=:circle, markersize=3.5, 
            color=my_colors,
            xlabel="Time [hr]", ylabel="Δ Height [μm/hr]",
            linewidth=2)
vline!([7.5], color=ColorSchemes.okabe_ito[1], linestyle=:dash, label=false, linewidth=3)
annotate!(5.5, 1, text("I", 12, "courier"))
annotate!(11, 1, text("II", 12, "courier"))
plot!(size=(400, 250))
#savefig("figs/fig2/fig2_time.pdf")
##
nr_profiles = npzread("data/timelapses/columns/profiles_nr.npy")
nr_heights = [NaNMath.mean(nr_profiles[i,:4000:6000]) for i=1:36]
nr_heights = reshape(nr_heights, (3,12))
r_profiles = npzread("data/timelapses/columns/profiles_r.npy")
r_heights = [NaNMath.mean(r_profiles[i,:2500:4500]) for i=1:36]
r_heights = reshape(r_heights, (3,12))
p3 = plot(nr_heights, color=1, marker=:circle)
plot!(r_heights, color=2, marker=:diamond)
plot!(legend=false, grid=false)
nr_h, nr_e = mean(nr_heights, dims=2), std(nr_heights, dims=2)
r_h, r_e = mean(r_heights, dims=2), std(r_heights, dims=2)
p3 = plot([nr_h, r_h], yerror=[nr_e r_e], marker=[:circle :diamond], markersize=5,  linewidth=3,color=[:gray ColorSchemes.okabe_ito[1]],
          label=["NR" "R"])
plot!(legend=(0.15, 0.25), ylim=(0, 320), xticks=[1,2,3], xlabel="Iteration", ylabel="Height [μm]", grid=false)
annotate!(1.15, 300, "C")

total_nr, total_r = maximum(nr_heights, dims=1), sum(r_heights, dims=1)
#=
bar!([mean(total_nr), mean(total_r)], yerror=[std(total_nr), std(total_r)], label=false,
     inset = (1, bbox(0.15, 0.6, 0.15, 0.35)),subplot=2, title="Tot. growth [μm]", color=[:gray, ColorSchemes.okabe_ito[1]],
     titlefontsize=8, xticks=[], yticks=[0, 250, 500])
=#
p4 = bar([mean(total_nr), mean(total_r)], yerror=[std(total_nr), std(total_r)], 
         label=false, xticks=([1,2], ["NR", "R"]), color=[:gray, ColorSchemes.okabe_ito[1]], 
         yticks=[], ylabel="Total Height [μm]", grid=false, ylim=(0, 550))
annotate!(0.94, 505, "D")
annotate!(1, 15, text("295.4 μm", 6, "courier", :left, rotation=90, color=:black))
annotate!(2, 15, text("462.1 μm", 6, "courier", :left, rotation=90, color=:black))

p5 = plot(x, myc, xlabel="Distance from interface", grid=false, color=:gray, linewidth=3, label="Concentration", size=(400, 400), xlim=(0,2))
plot!(x, sol_ana ./ maximum(sol_ana), color=ColorSchemes.okabe_ito[1], linewidth=3, label="Cumulative", xlim=(0,2))
plot!(x, sol_dis ./ maximum(sol_dis), color=ColorSchemes.okabe_ito[1],  linewidth=3,linestyle=:dash,label="Approximation")
vline!([sqrt(D*dt)], color=:black, alpha=0.6, linewidth=1.5, style=:dash, legend=:right, label=false, )
xticks!([0.5, 1.0, 2.0, 3.0]*sqrt(D*dt), ["0.5L", "L", "2L", "3L"], ylabel="Lim. nutrient (a.u.)", size=(400, 200), yticks=[])
annotate!(0.12, 0.96, "E")

#savefig("figs/fig2/2_conconly.svg")
l = @layout [[a{0.6h}
              b{0.7w} c{0.3w}] [d{0.6h}
                        e{0.4h}]]
plot(p1, p3, p4, p2, p5, layout=l, size=(700, 450), dpi=300)
savefig("figs/fig2/fig2.svg")
##
nr_bounds = npzread("data/timelapses/columns/bounds_A.npy")
nr_radius = reshape(nr_bounds[:,2]-nr_bounds[:,1], (3,12))*0.865
r_bounds = npzread("data/timelapses/columns/bounds_B.npy")
r_radius = reshape(r_bounds[:,2]-r_bounds[:,1], (3, 12))*0.865
plot(nr_radius, color=1)
plot!(r_radius, color=2, legend=false)
#
nr_volume = pi*(nr_heights .* nr_radius)
r_volume = pi*(r_heights .* r_radius)
plot(nr_volume, color=1)
plot!(r_volume, color=2)

nr_h, nr_e = mean(nr_volume, dims=2), std(nr_volume, dims=2)
r_h, r_e = mean(r_volume, dims=2), std(r_volume, dims=2)
p_v1 = plot([nr_h, r_h], yerror=[nr_e r_e], marker=[:circle :diamond], markersize=5,  
            linewidth=3,color=[:gray ColorSchemes.okabe_ito[1]],
            ylim=(0, 8e6),
            label=["NR" "R"])
plot!(legend=(0.15, 0.25), xticks=[1,2,3], xlabel="Iteration", ylabel="Volume [μm³]")
#
total_nr, total_r = maximum(nr_volume, dims=1), sum(r_volume, dims=1)
#=
bar!([mean(total_nr), mean(total_r)], yerror=[std(total_nr), std(total_r)], label=false,
     inset = (1, bbox(0.15, 0.6, 0.15, 0.35)),subplot=2, title="Tot. growth [μm]", color=[:gray, ColorSchemes.okabe_ito[1]],
     titlefontsize=8, xticks=[], yticks=[0, 250, 500])
=#
p_v2 = bar([mean(total_nr), mean(total_r)], yerror=[std(total_nr), std(total_r)], 
         label=false, xticks=([1,2], ["NR", "R"]), color=[:gray, ColorSchemes.okabe_ito[1]], 
         yticks=[], ylabel="Total Volume [μm³]", ylim=(0, 8.2e6))
annotate!(1, 5e5, text("6.82e6 μm³", 8, "courier", :left, rotation=90, color=:black))
annotate!(2, 5e5, text("7.11e6 μm³", 8, "courier", :left, rotation=90, color=:black))
l = @layout [a{0.7w} b]
plot(p_v1, p_v2, grid=false, size=(500, 250), layout=l, bottom_margin=3mm, left_margin=3mm, dpi=500)
savefig("figs/fig2/fig2_volumes.svg")

##
mean(total_r)