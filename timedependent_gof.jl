#=
Here we explore how the fits change depending on the amount of data. 
The main subject of study is JT305, since we managed to obtain 88 hours
of timelapse data across 3 replicates. 

From times 40h to 88h, we calculate the best fits and also 1000
bootstrapped parameters to obtain Confidence intervals (CI).

The data suggests that β is the last parameter to 'narrow down'
and therefore h_max will be imprecise until this converges.
=#
using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Plots, StatsPlots 

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function fit_data(t_data, h_data, model=interface, pguess=[0.8, 0.1, 15.0])
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
    result_ode = DiffEqFlux.sciml_train(loss, pguess, maxiters=500,
                                        lower_bounds = [1e-3, 1e-4, 1],
                                        upper_bounds = [1e1, 1, 1e2])
    return result_ode
end 
function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size-1) for x in sample_start])
    si = [(x % size(df)[1])+1 for x in sample_temp]
    return df[si,:]  
end
function boot_fit(df, n)
    p_bootstrap = []
    n_boots = Int(floor(length(df.time)/(3*5)))
    while length(p_bootstrap) < n
        Threads.@threads for i=1:n-length(p_bootstrap)
            boot_df = block_bootstrap(df, n_boots, 5)
            try
                myfit = fit_data(boot_df.time, boot_df.avg_height)
                append!(p_bootstrap, [myfit.u])
            catch
            end
        end
    end
    return reduce(hcat, p_bootstrap)'
end
df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== "jt305" && x.replicate in ["A", "B", "C"], df)
time_reference = filter(x-> x.replicate in ["C"], df).time
##
my_df = DataFrame("α"=>Float32[], "β"=>Float32[], "L"=>Float32[], "id"=>String[],
                  "time"=>Float32[], "h_sample"=>Float32[], "strain"=>String[])
n = 1000
for i=74:128
    println(i)
    tf = filter(x-> x.time .<= time_reference[i], df)
    best_fit = fit_data(tf.time, tf.avg_height).u
    boot_params = boot_fit(tf, n)
    id_names = lpad.(Array(1:n),3,"0")
    boot_data = vcat([best_fit ... "best"], hcat(boot_params, id_names))
    time_list = repeat([maximum(tf.time)], n+1)
    height_list = repeat([maximum(tf.avg_height)], n+1)
    strain_list = repeat(["jt305"], n+1)
    boot_data = [boot_data time_list height_list strain_list]
    [push!(my_df, boot_data[i,:]) for i=1:size(boot_data)[1]]
end
##
CSV.write("data/sims/bootstrap/time_jt305.csv", my_df)
