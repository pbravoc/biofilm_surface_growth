#= This code loops over the database, and returns the best fits
for the 48 hour data. One for each timelapse, and one for the aggregated
data.
It also calculated the best fit using 48h data + measurements from 
longtime_data.csv
=#
using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob
##
function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size) for x in sample_start])
    si = [(x % size(df)[1])+1 for x in sample_temp]
    return df[si,:]  
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

function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function fit_data(t_data, h_data, model) 
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    prob = ODEProblem(model["equation"], model["u0"], (0.0, t_data[end]), model["p_guess"])
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Only save height
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, model["p_guess"], maxiters=500,
                                        lower_bounds = model["p_low"],
                                        upper_bounds = model["p_high"])
    return result_ode
end 

## Define the models
interface_model = Dict("name"=>"interface",
                       "equation"=>interface, 
                       "u0"=>[0.1],
                       "p_guess"=>[0.8, 0.1, 15.0],
                       "p_low"=>[1e-5, 1e-5, 1e0],
                       "p_high"=>[1e3, 1e2, 1e5])

logisticnd_model = Dict("name"=>"logisticnd",
                        "equation"=>logistic_n,
                        "u0"=>[0.1, 1.0],
                        "p_guess"=>[0.1, 300.0, 0.5, 0.1],
                        "p_low"=>[1e-5, 1e0, 1e-5, 1e-5],
                        "p_high"=>[1e5, 1e5, 1e0, 1e5])

logistic_model = Dict("name"=>"logistic",
                      "equation"=>logistic, 
                      "u0"=>[0.1],
                      "p_guess"=>[0.1, 300.0],
                      "p_low"=>[1e-5, 1e0],
                      "p_high"=>[1e3, 1e5])

nutrient_model = Dict("name"=>"nutrient",
                      "equation"=>nutrient_n, 
                      "u0"=>[0.1, 1.0],
                      "p_guess"=>[1.0, 1e-1, 0.5, 0.1],
                      "p_low"=>[1e-5, 1e-5, 1e-5, 1e-5],
                      "p_high"=>[1e3, 1e2, 1e3, 1e5])
##
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0 && x.time .<48 && x.replicate in ["A", "B", "C"], Df)                    # Smaller than 0 values don't make physical sense
df2 = DataFrame(CSV.File("data/timelapses/longtime_data.csv"))
strain_list = unique(Df.strain)
P = []
Strain = []
Fit = []
model_choice = nutrient_model
# Get the best fits for less than 48h 
for strain in strain_list 
    println(strain)
    df = filter(x-> x.strain .== strain && x.time .<48, Df)              
    fit_params = fit_data(df.time, df.avg_height, model_choice)
    append!(P, [fit_params.u])
    append!(Strain, [strain])
    append!(Fit, ["48h"])
end

##
pf = hcat(DataFrame("strain"=>Strain, "fit"=>Fit),#, "u0"=>U),
          DataFrame(Matrix(reduce(hcat, P)'), :auto))
# Save to file
CSV.write("data/timelapses/parameters/params_"*model_choice["name"]*".csv", pf)

## If you want more detailed fits (across strains and long times)
#=


##
# Get the best fits for each timelapse
for strain in strain_list 
    println(strain)
    for replicate in ["A", "B", "C"]     
        df = filter(x-> x.strain .== strain && x.time .<48 &&
                        x.replicate .==replicate, Df)
        fit_params = fit_data(df.time, df.avg_height, model_choice)
        append!(P, [fit_params.u])
        append!(Strain, [strain])
        append!(Fit, [replicate])
    end
end
# This is to get the 'long time' best fit.
strain_list = ["bgt127", "gob33", "jt305"]
df = filter(x-> x.strain in strain_list && x.time .< 48 &&
                x.replicate in ["A", "B", "C"], Df)
for strain in strain_list 
    tf = filter(x-> x.strain .== strain, df)
    tf2 = filter(x-> x.strain .== strain, df2)
    t = append!(tf.time, tf2.time)
    h = append!(tf.avg_height, tf2.avg_height)
    p = fit_data(t, h, model_choice)
    append!(P, [p.u])
    append!(Strain, [strain])
    append!(Fit, ["long"])
end
=#