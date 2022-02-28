using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes
using DifferentialEquations, DiffEqFlux

function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    return t, h_avg, h_std
end
G(z, zstar) = z < zstar ? z : zstar 
function nutrient_n(du, u , p, t)
    h, c = u
    α, β, K_c, ϵ = p         
    du[1] = α*h*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
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
function param_fit(model_ode, model_params, pl, ph, u0)
    prob = ODEProblem(model_ode, u0, (0.0, 50.0), model_params) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_avg)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, params_guess[i], 
                                        lower_bounds = pl,
                                        upper_bounds = ph, 
                                        maxiters=100)
    return result_ode 
end

## Unbounded parameters fitting and simulation
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.time .<= 48, Df)
t, h_avg, h_std = get_average(df, "bgt127")
h_avg = abs.(h_avg)

u01, u02 = [0.2], [0.2, 1.0]
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
                [0.0, 0.0, 0.0]]
params_high =   [[Inf, 1e1, 1.0, 1.0], 
                [Inf, 1e4, 1.0, 1.0],
                [Inf, 1e1, 5e2, 1.0, 1.0],
                [Inf, 1e4],
                [30.0, 1e1, 50.0]]
model_params = []
solutions = []
sol_guess = []
for i=1:5
    print(models[i])
    problem = ODEProblem(models[i], starting_conditions[i], 
                         (0.0,48.0), params_guess[i])
    sol = solve(problem, saveat=0.5, save_idxs=1)
    append!(sol_guess, [sol])
end
plot()
[plot!(s) for s in sol_guess]
scatter!(t, h_avg)
   
##
for i=1:length(models)
    print(models[i])
    p = param_fit(models[i], params_guess[i], params_low[i],
                  params_high[i], starting_conditions[i])
    problem = ODEProblem(models[i], starting_conditions[i], (0.0,48.0), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    append!(model_params, [p.u])
    append!(solutions, [sol])
end

##
model_params 
##
for i=1:5
    prob = ODEProblem(models[i], starting_conditions[i], (0.0, 50.0), params_guess[i]) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_avg)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, params_guess[i], 
                                        lower_bounds = params_low[i],
                                        upper_bounds = params_high[i], 
                                        maxiters=100)

    probnew = ODEProblem(models[i], starting_conditions[i], (0.0, 50.0), result_ode)
    sol = solve(probnew, saveat=t, save_idxs=1)
    append!(solutions, [sol])
    append!(model_params, [result_ode])
end
    #scatter(t, h_avg)
    #plot!(sol)

##
pf = DataFrame()
pf.names = ["Nutrient_n", "Logistic_n", 
            "Interface_n", "Nutrient", "Interface"]
pf.parameters = [x.u for x in model_params] 
CSV.write("data/sims/bgt127_params.csv", pf)
##
##
pf = DataFrame("time"=>t, "data"=>h_avg, "data_error"=>h_std,
               "nutrient_n"=>solutions[1].u, "logistic_n"=>solutions[2].u, 
               "interface_n"=>solutions[3].u, "logistic"=>solutions[4].u, 
               "interface"=>solutions[5].u)
CSV.write("data/sims/bgt127.csv", pf)

##
pf = DataFrame("time"=>t, "data"=>h_avg, "data_error"=>h_std,
               "nutrient_n"=>solutions[1].u, "logistic_n"=>solutions[2].u, 
               "interface_n"=>solutions[3].u, "logistic"=>solutions[4].u, 
               "interface"=>solutions[5].u)
CSV.write("data/sims/bgt127.csv", pf)