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
df = filter(x-> x.replicate in ["A", "B", "C"] && x.strain .== "jt305", Df)[1:end-1]
@df df plot(:time, :avg_height, group=(:replicate))
##
t, h_avg, h_std = get_average(df, "jt305")
h_avg = abs.(h_avg)
##
u01, u02 = [0.2], [0.2, 1.0]
models = [nutrient_n, logistic_n, interface_n, logistic, interface]
starting_conditions = [u02, u02, u02, u01, u01]
params_guess = [[0.8, 0.001, 0.9, 0.12], 
                [0.5, 100, 0.2, 0.05],
                [0.8, 0.07, 15, 0.01, 0.01],
                [0.5, 150.0],
                [0.8, 0.07, 15]]
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
model_params = []
solutions = []
##
for i=1:5
    print(models[i])
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
CSV.write("data/sims/gob33_params.csv", pf)
##
##
pf = DataFrame("time"=>t, "data"=>h_avg, "data_error"=>h_std,
               "nutrient_n"=>solutions[1].u, "logistic_n"=>solutions[2].u, 
               "interface_n"=>solutions[3].u, "logistic"=>solutions[4].u, 
               "interface"=>solutions[5].u)
CSV.write("data/sims/gob33.csv", pf)

##
@df pf scatter(:time, :data, yerror=:data_error, color=:black, alpha=0.3)
@df pf plot!(:time, :interface, linestyle=:dash, color=1, linewidth=2)
@df pf plot!(:time, :interface_n, color=1, linewidth=2)
@df pf plot!(:time, :logistic_n, color=2, linewidth=2)
@df pf plot!(:time, :logistic, linestyle=:dash, color=2, linewidth=2)
@df pf plot!(:time, :nutrient_n, color=3, linewidth=2)


##
strain_name = "gob33"
parameters = DataFrame(CSV.File("data/sims/"*strain_name*"_params.csv"))
data = DataFrame(CSV.File("data/sims/"*strain_name*".csv"))
parameters = [parse.(Float64, split(parameters.parameters[i][2:end-1], ", ")) for i=1:5]
solutions = []
for i=4:5
    prob = ODEProblem(models[i], starting_conditions[i], (0.0, 350.0), parameters[i]) # Set the problem
    sol = solve(prob, saveat=0.5, save_idxs=1)
    append!(solutions, [sol.u])
end
plot(sol.t, solutions)
pf = DataFrame("time"=>sol.t, "logistic"=>solutions[1], "interface"=>solutions[2])
CSV.write("data/sims/"*strain_name*"_lt.csv", pf)
#models = [logistic, interface]

##
plot(sol.t, solutions)
