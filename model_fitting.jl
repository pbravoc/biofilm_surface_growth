using DifferentialEquations, DiffEqFlux, Plots
using DataFrames, CSV 
using StatsPlots

## Models
G(z, zstar) = z < zstar ? z : zstar 
function nutrient_n(du, u , p, t)
    h, c = u
    α, β, K_c, ϵ = p         
    du[1] = α*h*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end
function nutrient(du, u , p, t)
    h = u[1]
    α, β = p         
    du[1] = α*h - β*h 
    return du 
end
function logistic_n(du, u , p, t)
    h, c = u
    α, K_h, K_c, ϵ = p         
    du[1] = α*h*(1- h/(K_h))*(c/(K_c + c)) 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end
function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    #du[1] = α * h - (α*h*h/K_h + h)
    du[1] = α*h*(1- h/(K_h))
    #du[1] = α * h *(1- (h/(K_h + h)))
    return du 
end
function interface_n(du, u, p, t)
    h, c = u 
    α, β, hstar, K_c, ϵ = p
    du[1] = α*G.(h, hstar)*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du
end
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function param_fit(model_ode, model_params, u0)
    prob = ODEProblem(model_ode, u0, (0.0, 50.0), model_params) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=df.time, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- df.loess_height)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, model_params,
                                        maxiters=100)
    return result_ode 
end
## Nutrient depletion + n 
u0 = [0.3, 1.0]
p = [0.8, 0.001, 0.9, 0.12]
prob = ODEProblem(nutrient_n, u0, (0.0, 50.0), p)
sol = solve(prob, save_idxs=1)
plot(sol)

## Nutrient depletion
u0 = [0.1]
p = [0.6*0.666, 0.1]
prob = ODEProblem(nutrient, u0, (0.0, 50.0), p)
sol = solve(prob)
plot!(sol)

## Logistic_n 
u0 = [0.3, 1.0]
p = [0.7, 200, 0.2, 0.05]
prob = ODEProblem(logistic_n, u0, (0.0, 50.0), p)
sol = solve(prob)
plot(sol)

## logistic
u0 = [0.1]
p = [0.7, 200.0]
prob = ODEProblem(logistic, u0, (0.0, 50.0), p)
sol = solve(prob)
plot!(sol)

## Inferface limited + nutrients 
u0 = [0.1, 1.0]
p = [0.9, 0.07, 15, 0.01, 0.01]
prob = ODEProblem(interface_n, u0, (0.0, 50.0), p)
sol = solve(prob)
plot(sol)

## Inferface limited
u0 = [0.1]
p = [0.9, 0.07, 15]
prob = ODEProblem(interface, u0, (0.0, 50.0), p)
sol = solve(prob)
plot(sol)

## Loading data
my_strain, my_replicate = "bgt127", "A"
df =  filter(row -> row.replicate .== my_replicate && 
             row.strain .== my_strain, 
             DataFrame(CSV.File("data/timelapses/database_updated.csv")));
@df df scatter(:time, :loess_height, label=false)
 
##
p_nutrients_n = param_fit(nutrient_n, [0.8, 0.001, 0.9, 0.12], [df.loess_height[1], 1.0])
p_logistic_n = param_fit(logistic_n, [0.7, 200, 0.2, 0.05], [df.loess_height[1], 1.0])
p_interface_n = param_fit(interface_n, [0.9, 0.07, 15, 0.01, 0.01], [df.loess_height[1], 1.0])
p_logistic = param_fit(logistic, [0.7, 200.0], [df.loess_height[1]])
p_interface = param_fit(interface, [0.9, 0.07, 15], [df.loess_height[1]])
##
u01 = [df.loess_height[1]]
u02 = [df.loess_height[1], 1.0]
problem_nutrient_n = ODEProblem(nutrient_n, u02, (0.0, 350.0), p_nutrients_n)
sol_nutrient_n = solve(problem_nutrient_n, saveat=0.5, save_idxs=1)
problem_logistic_n = ODEProblem(logistic_n, u02, (0.0, 350.0), p_logistic_n)
sol_logistic_n = solve(problem_logistic_n, saveat=0.5, save_idxs=1)
problem_interface_n = ODEProblem(interface_n, u02, (0.0, 350.0), p_interface_n)
sol_interface_n = solve(problem_interface_n, saveat=0.5, save_idxs=1)
problem_logistic = ODEProblem(logistic, u01, (0.0, 350.0), p_logistic)
sol_logistic = solve(problem_logistic, saveat=0.5)
problem_interface = ODEProblem(interface, u01, (0.0, 350.0), p_interface)
sol_interface = solve(problem_interface, saveat=0.5)
##
@df df scatter(:time, :avg_height, color=:black, alpha=0.25, label=false)
#@df Df scatter(:time, :avg_height, color=:black, alpha=0.25, label=false)
plot!(sol_nutrient_n, color=1, linewidth=1.5, label="Birth-death")
plot!(sol_logistic_n, color=2,linewidth=1.5, label="Logistic")
plot!(sol_interface_n, color=3, linewidth=1.5,label="Interface",legend=:bottomright)
#plot!(sol_logistic, color=2,linewidth=1.5, linestyle=:dash, label="Logistic(nn)")
#plot!(sol_interface, color=3, linewidth=1.5,linestyle=:dash, 
      #label="Interfae(nn)", legend=:bottomright)
plot!(xlim=(0, 50), ylim=(-1, 220), xlabel="Time [hr]", ylabel="Height [μm]")
savefig("figs/timelapses/bgt127/fit_models_n.svg")
##
Df =  filter(row ->  row.strain .== my_strain && row.time .> 300,
             DataFrame(CSV.File("data/timelapses/database_updated.csv")));
@df Df scatter(:time, :avg_height, color=:black, alpha=0.25, label=false)
