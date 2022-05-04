#= 
Collection of all the useful functions for the project
=#
using DifferentialEquations, DiffEqFlux
using Plots, StatsPlots, ColorSchemes, Plots.PlotMeasures
using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using Glob

"""Re-samples your dataset in blocks"""
function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size) for x in sample_start])
    si = [(x % size(df)[1])+1 for x in sample_temp]
    return df[si,:]  
end

"""Logistic growth with nutrients"""
function ode_logistic_n(du, u , p, t)
    h, c = u
    α, K_h, K_c, ϵ = p         
    du[1] = α*h*(1- h/(K_h))*(c/(K_c + c)) 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end

"""Logistic growth without nutrients"""
function ode_logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    du[1] = α*h*(1- h/(K_h))
    return du 
end

"""Structure that holds a model with initial 
conditions, boundaries and parameter guesses"""
struct model_struct
    u0::Array
    p_guess::Array
    p_low::Array 
    p_high::Array 
    ode::Function
    name::String
end

function fit_data(t_data, h_data, mstruct, bounded=false)
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    u0 = mstruct.u0
    prob = ODEProblem(mstruct.ode, u0, (0.0, t_data[end]), mstruct.p_guess)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    if bounded 
        result_ode = DiffEqFlux.sciml_train(loss, mstruct.p_guess, 
                                        lower_bounds = mstruct.p_low,
                                        upper_bounds = mstruct.p_high, 
                                        maxiters=100)
    else
        result_ode = DiffEqFlux.sciml_train(loss, mstruct.p_guess)
    end
    return result_ode
end 