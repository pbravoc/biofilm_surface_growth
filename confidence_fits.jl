using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using Plots, StatsPlots
using DifferentialEquations, DiffEqFlux

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u 
    α, β, hstar = p[:,1], p[:,2], p[:,3]
    for i in 1:length(u)
        du[i] = α[i] * G(h[i], hstar[i]) - β[i] * h[i] 
    end
end

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))

function plot_longtime(strain_name, c)
    df = filter(x-> x.strain .== strain_name && x.replicate in unique(Df.replicate)[4:end], Df)
    pf = DataFrame(CSV.File("data/sims/bootstrap/boot_"*strain_name*".csv"))
    p = Array(pf[1:1000, 1:3])
    u0 = repeat([0.2], 1000)
    prob = ODEProblem(interface, u0, (0.0, 348.0), p)
    sol = solve(prob, saveat=1.0)
    quant_heights = reduce(hcat, [quantile(sol[i], [0.025, 0.25, 0.75, 0.975]) for i=1:length(sol.t)])
    mean_height = [median(sol[i]) for i=1:length(sol.t)]
    plot!(sol.t, mean_height, linewidth=2, color=c, label=false)
    #plot!(quant_heights[2,:], alpha=0.0, fillrange = quant_heights[3,:], fillalpha = 0.15, c = c, label = false)
    plot!(quant_heights[1,:], alpha=0.0, fillrange = quant_heights[4,:], fillalpha = 0.15, c = c, label = false,)
    @df df scatter!(:time, :avg_height, color=c, legend=:topleft, label=strain_name)
end
plot(xlabel="Time [hr]", ylabel="Height [μm]", dpi=500)
plot_longtime("bgt127", 1)
plot_longtime("jt305", 2)
plot_longtime("gob33", 3)
savefig("figs/bootstrapping/longtime_95CI.png")
