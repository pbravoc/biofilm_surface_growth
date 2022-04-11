using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob

function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size) for x in sample_start])
    si = [(x % size(df)[1])+1 for x in sample_temp]
    return df[si,:]  
end
G(z, zstar) = z < zstar ? z : zstar 

function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

function fit_data(t_data, h_data, model, pguess=[0.8, 100.0])
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    u0 = [0.1]
    prob = ODEProblem(model, u0, (0.0, t_data[end]), pguess)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, pguess)
    return result_ode
end 

function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    du[1] = α*h*(1- h/(K_h))
    return du 
end

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0, Df)                    # Smaller than 0 values don't make physical sense
strain_list = unique(Df.strain)
model_choice, n_parameters = logistic, 2
P = []
Strain = []
Fit = []
##
# Get the best fits for all timepoints
for strain in strain_list 
    df = filter(x-> x.strain .== strain, Df)              
    fit_params = fit_data(df.time, df.avg_height, model_choice)
    append!(P, [fit_params.u])
    append!(Strain, [strain])
    append!(Fit, ["all"])
end
# Get the best fits for less than 48h 
for strain in strain_list 
    df = filter(x-> x.strain .== strain && x.time .<48, Df)              
    fit_params = fit_data(df.time, df.avg_height, model_choice)
    append!(P, [fit_params.u])
    append!(Strain, [strain])
    append!(Fit, ["48h"])
end
# Get the best fits for each timelapse
for strain in strain_list 
    for replicate in ["A", "B", "C"]     
        df = filter(x-> x.strain .== strain && x.time .<48 &&
                        x.replicate .==replicate, Df)
        fit_params = fit_data(df.time, df.avg_height, model_choice)
        append!(P, [fit_params.u])
        append!(Strain, [strain])
        append!(Fit, [replicate])
    end
end

pf = hcat(DataFrame("strain"=>Strain, "fit"=>Fit),
          DataFrame(Matrix(reduce(hcat, P)'), :auto))
## Save to file
#parameter_df.h_max = parameter_df.α .* parameter_df.L ./ parameter_df.β
CSV.write("data/timelapses/fit_params_"*string(model_choice)*".csv", pf)