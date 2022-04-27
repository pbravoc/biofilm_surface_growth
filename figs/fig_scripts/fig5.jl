using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors
using LsqFit
using DifferentialEquations

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

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

function interface_sim(p)
    prob = ODEProblem(interface, [0.1], (0.0, 48.0), p)
    sol = solve(prob, saveat=0.5, save_idxs=1)
    return sol
end

function order_average(Df, strain_name, subplot_id)
    df = filter(x->x.strain .== strain_name,Df)
    n = maximum(df.order)
    t, t_e, h, h_e = zeros(n), zeros(n), zeros(n), zeros(n)
    for i=1:maximum(df.order)
        tf = filter(x->x.order .== i,df)
        t[i], t_e[i] = mean(tf.time), std(tf.time)
        h[i], h_e[i] = mean(tf.avg_height), std(tf.avg_height)
    end
    scatter!(t, h, xerror=t_e, yerror=h_e, color=:black, markersize=1, alpha = 0.5,subplot=subplot_id)
end
##
df = DataFrame(CSV.File("data/timelapses/database.csv"))
pf = DataFrame(CSV.File("data/timelapses/fit_params.csv"))
order = []
for strain in unique(df.strain)
    for repli in unique(df.replicate)
        tf = filter(row->row.replicate.==repli && row.strain .== strain, df);
        my_order = Array(1:size(tf)[1])
        append!(order, my_order)  
    end  
end
df.order = order
df = filter(x-> x.time .< 48 && x.avg_height .> 0 && 
                x.replicate in ["A", "B", "C"], df)  
## 
solutions = []
for i=1:size(pf)[1]
    sol = interface_sim(pf[i,3:5])
    append!(solutions, [sol])
end
pf.solutions = solutions
##
my_plots = []
strain_names = unique(df.strain)
P = plot(layout=(2,3), size=(800, 500))
for i=1:6
    strain = strain_names[i]
    order_average(df, strain, i)
    best_sol = filter(x-> x.strain .== strain && 
                x.fit .== "48h", pf)[1,7]
    plot!(best_sol, linewidth=2, subplot=i, color=:red)
    plot!(xlabel="Time [hr]", ylabel="Height [μm]", title=strain, subplot=i,
    legend=false)
end 
display(P)
savefig("figs/fig4/fig4.svg")
##
p1 = @df filter(x-> x.fit .== "48h", pf) scatter(:α, :β, group=:strain, markersize=7)
@df filter(x-> x.fit .== "all", pf) scatter!(:α, :β, group=:strain, marker=:diamond,label=false, c=[1 2 3 4 5 6], markersize=7)
@df filter(x-> x.fit in ["A", "B", "C"], pf) scatter!(:α, :β, group=:strain, label=false, c=[1 2 3 4 5 6], alpha=0.5)
plot!(xlabel="α [μm/hr]", ylabel="β [μm/hr]", legend=false)

p2 = @df filter(x-> x.fit .== "48h", pf) scatter(:β, :L, group=:strain, markersize=7)
@df filter(x-> x.fit .== "all", pf) scatter!(:β, :L, group=:strain, marker=:diamond,label=false, c=[1 2 3 4 5 6], markersize=7)
@df filter(x-> x.fit in ["A", "B", "C"], pf) scatter!(:β, :L, group=:strain, label=false, c=[1 2 3 4 5 6], alpha=0.5)
plot!(xlabel="β [μm/hr]", ylabel="L [μm]", legend=false)

p3 = @df filter(x-> x.fit .== "48h", pf) scatter(:α, :L, group=:strain, markersize=7)
@df filter(x-> x.fit .== "all", pf) scatter!(:α, :L, group=:strain, marker=:diamond,label=false, c=[1 2 3 4 5 6], markersize=7)
@df filter(x-> x.fit in ["A", "B", "C"], pf) scatter!(:α, :L, group=:strain, label=false, c=[1 2 3 4 5 6], alpha=0.5)
plot!(xlabel="α [μm/hr]", ylabel="L [μm]")
scatter!([2.0], [15.0], label="48h", color=:black)
scatter!([2.0], [15.0], label="All", color=:black, marker=:diamond)
scatter!([2.0], [15.0], label="Replicate", color=:black, alpha=0.5)
plot!(xlim=(0.27, 1.5))
plot(p1, p2, p3, layout=(1,3), size=(900, 300), left_margin=4mm, bottom_margin=5mm)
savefig("figs/fig4/best_fits.svg")

##
p1 = @df filter(x-> x.fit .== "48h", pf) scatter(:α, :β, :L, group=:strain, markersize=5)
@df filter(x-> x.fit .== "all", pf) scatter!(:α, :β, :L, group=:strain, marker=:diamond,label=false, c=[1 2 3 4 5 6], markersize=5)
@df filter(x-> x.fit in ["A", "B", "C"], pf) scatter!(:α, :β, :L, group=:strain, label=false, c=[1 2 3 4 5 6], alpha=0.8, markersize=3)
plot!(xlabel="Alpha [um/hr]", ylabel="Beta [um/hr]", zlabel="L [um]", size=(500, 400))
savefig("figs/fig4/best_fits.html")
