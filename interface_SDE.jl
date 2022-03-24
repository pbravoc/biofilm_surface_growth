using DataFrames, CSV
using DifferentialEquations
using Plots, StatsPlots
using SpecialFunctions

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function σ_interface(du, u , p ,t)
    h = u[1] 
    α, β, hstar = p
    du[1] = 0.1*α*G.(h, hstar) + 0.1*β*h
    return du
end

function sim_trajectories(p, s)
    u0 = [0.1]
    prob_sde = SDEProblem(interface,σ_interface, u0, (0.0, 48.0), p)
    ensembleprob = EnsembleProblem(prob_sde)
    sol = solve(ensembleprob,EnsembleThreads(),trajectories=3, saveat=0.5, save_idxs=1)
    all_t = reduce(vcat, [x.t for x in sol])
    all_h = reduce(vcat, [x.u for x in sol])
    α = repeat([p[1]], length(all_t))
    β = repeat([p[2]], length(all_t))
    L = repeat([p[3]], length(all_t))
    r = reduce(vcat, [repeat(["A"], Int(length(all_t)/3)),
                      repeat(["B"], Int(length(all_t)/3)),
                      repeat(["C"], Int(length(all_t)/3))])
    group_name = repeat([s], length(all_t))
    df = DataFrame("S"=>group_name, "r"=>r, "α"=>α, "β"=>β, "L"=>L,
                   "t"=>all_t, "h"=>all_h)
    return df
end

p = [0.8, 0.1, 15.0]
u0 = [0.1]
prob_sde = SDEProblem(interface,σ_interface, u0, (0.0, 48.0), p)
ensembleprob = EnsembleProblem(prob_sde)
sol = solve(ensembleprob,EnsembleThreads(),trajectories=3, saveat=0.5, save_idxs=1)
scatter(sol)
##
α_range = 0.3:0.1:0.7
β_range = 0.05:0.01:0.09
L_range = 10.0:10.0:50.0

s = 0
Df = DataFrame("S"=>Int64[],"r"=>String[], "α"=>Float64[],
               "β"=>Float64[], "L"=>Float64[], "t"=>Float64[], 
               "h"=>Float64[]) 
for a in α_range, b in β_range, l in L_range
    p = [a, b, l]
    df = sim_trajectories(p, s)
    append!(Df, df)
    s += 1
    print(s)
end
##
@df Df scatter(:t, :h, group=(:S), label=false, alpha=0.3, markersize=2, markerstrokecolor=:auto)
##
CSV.write("data/sims/simulated_data.csv", Df)
##
@df filter(x-> x.S==100, Df) scatter(:t, :h)
##
@df pf scatter(:α, :β, :L, group=(:S), label=false, alpha=0.9, markersize=2, markerstrokecolor=:auto, xlim=(0.2, 1.0), ylim=(0.0, 0.2))

