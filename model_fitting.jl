using DifferentialEquations, DiffEqFlux, Plots
using DataFrames, CSV 
using StatsPlots
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
function param_fit_box(model_ode, model_params, p_l, p_h, u0)
    prob = ODEProblem(model_ode, u0, (0.0, 50.0), model_params) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=df.time, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- df.loess_height)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, model_params, lower_bounds=p_l, upper_bounds=p_h,
                                        maxiters=100)
    return result_ode 
end
function d_height(s)
    h_change = zeros(length(s.u)) 
    h_change .= NaN
    if length(s.u) > 2
        h_change[1:end-1] = (s.u[2:end]-s.u[1:end-1]) ./ 
                            (s.t[2:end]-s.t[1:end-1])
    end
    return h_change
end
my_strain, my_replicate = "bgt127", "A"
df =  filter(row -> row.replicate .== my_replicate && 
             row.strain .== my_strain, 
             DataFrame(CSV.File("data/timelapses/database_updated.csv")));
@df df scatter(:time, :loess_height, label=false)
 
##
u01, u02 = [df.loess_height[1]], [df.loess_height[1], 1.0]
models = [nutrient_n, logistic_n, interface_n, logistic, interface]
starting_conditions = [u02, u02, u02, u01, u01]
params_guess = [[0.8, 0.001, 0.9, 0.12], 
                [0.7, 200, 0.2, 0.05],
                [0.9, 0.07, 15, 0.01, 0.01],
                [0.7, 200.0],
                [0.9, 0.07, 15]]
params_low =    [[1e-3, 1e-4, 0.0, 0.0], 
                [1e-2, 10.0, 0.0, 0.0],
                [1e-2, 1e-5, 5.0, 0.0, 0.0],
                [1e-3, 50.0],
                [1e-3, 1e-3, 5.0]]
params_high =   [[1e3, 1e1, 1.0, 1.0], 
                [1e3, 1e3, 1.0, 1.0],
                [1e3, 1e1, 5e2, 1.0, 1.0],
                [1e3, 1e3],
                [1e3, 1e2, 5e2]]
model_params = []
solutions = []
for i=1:length(models)
    print(models[i])
    p = param_fit_box(models[i], params_guess[i], params_low[i], params_high[i], 
                      starting_conditions[i])
    print(p)
    problem = ODEProblem(models[i], starting_conditions[i], (0.0,50.0), p)
    sol = solve(problem, saveat=0.5, save_idxs=1)
    append!(model_params, [p.u])
    append!(solutions, [sol])
end
##
colors = [1,2,3,2,3]
line_style = [:solid, :solid, :solid, :dash, :dash]
plot()
for i in 1:length(solutions)
    plot!(solutions[i], label=models[i], linewidth=1.5, color=colors[i], linestyle=line_style[i])
end
@df df scatter!(:time, :avg_height, color=:black, alpha=0.25, label=false, legend=:topleft)

##

