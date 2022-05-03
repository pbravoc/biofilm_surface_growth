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
function fit_data(t_data, h_data, pguess=[0.8, 0.1, 15.0])
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    u0 = [0.1]
    prob = ODEProblem(interface, u0, (0.0, t_data[end]), pguess)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, pguess)
    return result_ode
end 

function boot_fit(df, n)
    p_bootstrap = []
    best_fit = fit_data(df.time, df.avg_height)
    while length(p_bootstrap) < n
        Threads.@threads for i=1:n-length(p_bootstrap)
            boot_df = block_bootstrap(df, 20, 5)
            try
                myfit = fit_data(boot_df.time, boot_df.avg_height, best_fit)
                append!(p_bootstrap, [myfit.u])
                print(length(p_bootstrap))
            catch
            end
        end
    end
    return reduce(hcat, p_bootstrap)'
end
## This is to do bootstrap fitting on 48h data
strain_name = "ea387"
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== strain_name && x.time .< 48 &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
data = boot_fit(df, 1000)
parameters_frame = DataFrame("α"=>data[:,1], "β"=>data[:,2], "L"=>data[:,3])
parameters_frame.h_max = parameters_frame.α .* parameters_frame.L ./ parameters_frame.β
CSV.write("data/sims/bootstrap/boot_"*strain_name*".csv", parameters_frame)


## And to aggregate all the different strains on one file
strain_names = ["bgt127", "jt305", "gob33", "y55", "bh1514", "ea387"]
Data = DataFrame()
for sname in strain_names
    df = DataFrame(CSV.File("data/sims/bootstrap/boot_"*sname*".csv"))
    df.strain = repeat([sname], 1000)
    append!(Data, df)
end
CSV.write("data/sims/bootstrap/all_bootstrap.csv", Data)

## This is to get the 'long time' best fit.
strain_list = ["bgt127", "gob33", "jt305"]
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain in strain_list && x.time .< 48 &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
df2 = DataFrame(CSV.File("data/timelapses/longtime_data.csv"))
params = []
for strain in strain_list 
    tf = filter(x-> x.strain .== strain, df)
    tf2 = filter(x-> x.strain .== strain, df2)
    t = append!(tf.time, tf2.time)
    h = append!(tf.avg_height, tf2.avg_height)
    p = fit_data(t, h)
    append!(params, [p])
    print(strain)
end