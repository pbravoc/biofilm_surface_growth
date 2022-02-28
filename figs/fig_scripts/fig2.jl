using DataFrames, CSV
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

default(palette = [cgrad(:Dark2_3)[x] for x in [0.0, 0.05, 0.15, 0.45, 0.5, 0.55, 0.9, 0.95, 1.0]])
default(palette = palette(:default))

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

l = @layout [[a{0.7h}
              b{0.3h}] [c{0.7h}
                        d{0.3h}]]

p1 = @df df plot(:time, :avg_height, group=:replicate,
                fillalpha=0.3, ylabel="Height [μm]", xlabel="Time [hr]",
                alpha=1.0,  grid=false, 
                legend=:topleft,
                marker=:circle, markersize=2, markerstrokecolor=:auto,
                linewidth=2, ribbon=:std_height)
@df df plot!(:time, :smooth_height, group=:replicate,
             alpha=0.8, legend=false, 
             grid=false,
             marker=:circle, markersize=1, markerstrokecolor=:auto,
             inset = (1, bbox(0.5, 0.45, 0.5, 0.45)),
             linewidth=1.5, subplot=2, yscale=:log10)

p2 = @df df plot(:avg_height, :slope, group=:replicate,
            fillalpha=0.3, alpha=0.8,
            ribbon=:slope_error, grid=false,
            marker=:circle, markersize=3, markerstrokecolor=:auto,
            xlabel="Height [μm]", ylabel="Δ Height [μm/hr]",
            linewidth=2)
p2 = vline!([27.5], color=:black, linestyle=:dash, legend=false)
p2 = annotate!(20, 1, text("I", 8, "courier"))
p2 = annotate!(40, 1, text("II", 8, "courier"))

p3 = plot(x, myc, xlabel="Distance from interface", grid=false, color=1, linewidth=2, label="Concentration", size=(400, 400), xlim=(0,2))
p3 = plot!(x, sol_ana ./ maximum(sol_ana), color=2, linewidth=2, label="Cumulative", xlim=(0,2))
p3 = plot!(x, sol_dis ./ maximum(sol_dis), color=2,  linewidth=2,linestyle=:dash,label="Approximation")
p3 = vline!([sqrt(D*dt)], color=:black, alpha=0.6, linewidth=1.5, style=:dash, legend=:right, label=false, ylabel="Lim. nutrient")
p3 = xticks!([0.5, 1.0, 2.0, 3.0]*sqrt(D*dt), ["0.5L", "L", "2L", "3L"])



cube_data = [161.07 102.99 108.80
             179.58 92.46 92.48
             182.84 96.82 127.31
             232.57 101.17 73.23
             183.57 83.39 81.21
             166.51 79.03 59.80
             161.07 89.92 64.51
             200.99 108.07 NaN
             248.18 82.66 NaN]
cube_sorted = [248.18 232.57 200.99 183.57 182.84 179.58 166.51 161.07 161.07
              108.07 102.99 101.17 96.82 92.46 89.92 83.39 82.66 79.03
              127.31 108.8 92.48 81.21 73.23 64.51 59.8 NaN NaN]

p4 = groupedbar(cube_sorted, color=1, label=false, grid=false, xticks=[1,2,3], ylim=(0, 280), xlabel="Iteration", ylabel="Max. Height [μm]")

plot(p1, p3, p2, p4,layout=l, size=(700, 450), dpi=300)
#savefig("figs/figs_temp/fig2.png")

