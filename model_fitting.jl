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
function fit_data(t_data, h_data, model, model_params, pl, ph, u0)
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    prob = ODEProblem(model, u0, (0.0, t_data[end]), model_params)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, model_params) 
                                        #lower_bounds = pl,
                                        #upper_bounds = ph)
    return result_ode
end 

## Bounded parameters fitting
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.replicate in ["A", "B", "C"] && x.strain .== "bgt127", Df)
@df df plot(:time, :avg_height, group=(:replicate))
##
df.avg_height = abs.(df.avg_height)
t, h_avg, h_std = get_average(df, "bgt127")
scatter(t, h_avg, yerror=h_std, color=:black, alpha=0.5)
##
u01, u02 = [0.2], [0.2, 1.0]
models = [nutrient_n, logistic_n, interface_n, logistic, interface]
starting_conditions = [u02, u02, u02, u01, u01]
params_guess = [[0.8, 0.1, 0.9, 0.12], 
                [0.5, 100, 0.2, 0.05],
                [0.8, 0.07, 15, 0.01, 0.01],
                [0.5, 150.0],
                [0.8, 0.07, 15]]
params_low =    [[0.1, 0.01, 0.01, 0.001], 
                [1e-2, 10.0, 0.001, 0.001],
                [1e-2, 1e-5, 5.0, 0.001, 0.001],
                [1e-3, 10.0],
                [0.01, 0.05, 3.0]]
params_high =   [[5.0, 1e1, 1.0, 1.0], 
                [5.0, 1000.0, 1.0, 1.0],
                [5.0, 2.0, 1000.0, 1.0, 1.0],
                [5.0, 1000.0],
                [5.0, 1.0, 50.0]]
##
model_params = []
solutions = []
##
i = 5
result_ode = fit_data(t, h_avg, models[i], 
                        params_guess[i], params_low[i], params_high[i],
                        starting_conditions[i])
##
problem = ODEProblem(models[i], starting_conditions[i], 
                        (0.0,48.5), result_ode)
sol = solve(problem, saveat=t, save_idxs=1)   
##
append!(solutions, [sol])
append!(model_params, [result_ode])
##
plot()
[plot!(s) for s in solutions]
plot!(xlabel="Time [hr]", ylabel="Height [μm]", legend=:topleft)
##
pf = DataFrame()
pf.names = ["Nutrient_n", "Logistic_n", 
            "Interface_n", "Nutrient", "Interface"]
pf.parameters = [x.u for x in model_params] 
CSV.write("data/sims/fig3_params_unbounded.csv", pf)
##
pf = DataFrame("time"=>t, "data"=>h_avg, "data_error"=>h_std,
               "nutrient_n"=>solutions[1].u, "logistic_n"=>solutions[2].u, 
               "interface_n"=>solutions[3].u, "logistic"=>solutions[4].u, 
               "interface"=>solutions[5].u)
CSV.write("data/sims/f3a_heights_unbounded.csv", pf)

##
pf = DataFrame("time"=>solutions2[1].t,
               "nutrient_n"=>solutions2[1].u, "logistic_n"=>solutions2[2].u, 
               "interface_n"=>solutions2[3].u, "logistic"=>solutions2[4].u, 
               "interface"=>solutions2[5].u)
CSV.write("data/sims/f3a_heights_unboundedlong.csv", pf)

##
plot()
[plot!(s) for s in solutions2]
plot!(xlabel="Time [hr]", ylabel="Height [μm]", legend=:topleft)

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
