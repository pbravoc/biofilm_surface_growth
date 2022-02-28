##
2+3
##
using DataFrames, CSV
using Statistics
using Plots, StatsPlots
using DifferentialEquations, DiffEqFlux

function logistic(du, u , p, t)
    h, c = u
    α, K_h, K_c, ϵ = p         
    du[1] = α*h*(1- h/(K_h))*(c/(K_c + c)) 
    du[2] = -ϵ*(h*c/(K_c + c)) 
    return du 
end

function nutrient_n(du, u , p, t)
    h, c = u
    α, β, K_c, ϵ = p         
    du[1] = α*h*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end

##
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && x.replicate =="A", Df)
p = [0.85, 0.001, 0.9, 0.12]
u0 = [0.2, 1.0]
@df df scatter(:time, :avg_height)

problem = ODEProblem(nutrient_n, u0, (0.0,48.0), p)
sol = solve(problem, saveat=0.5, save_idxs=1)
plot!(sol)

##
function loss(p)
    sol = solve(problem, Tsit5(), p=p, saveat=df.time, save_idxs=1) # Force time savings to match data
    sol_array = reduce(vcat, sol.u)
    loss = sum(abs2, sol_array .- df.avg_height)
    return loss, sol
end
result_ode = DiffEqFlux.sciml_train(loss, p,
                                    lower_bounds = [0.0, 1e-4, 0.0, 0.0], 
                                    upper_bounds = [Inf, Inf, Inf, 0.2],
                                    maxiters=100)

##
sol_fit = solve(ODEProblem(nutrient_n, u0, (0.0,548.0), result_ode), save_idxs=1)
plot(sol_fit)
@df df scatter!(:time, :avg_height)

##
function monod(du, u, p ,t)
    x, s = u 
    μ, K_s, ϵ = p 
    du[1] = μ*ϵ*x*s/(K_s + s)
    du[2] = -μ*s*(sum(u0)-s)/(K_s + s)
    return du 
end

p = [2.0, 3.0, 1.0]
u0 = [0.1, 1.0]
problem = ODEProblem(monod, u0, (0.0, 10.0), p)
sol = solve(problem, save_idxs=1)
plot(sol)