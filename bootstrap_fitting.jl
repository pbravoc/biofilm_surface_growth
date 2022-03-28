using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux

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
function fit_data(t_data, h_data)
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    u0 = [0.1]
    pguess = [0.8, 0.1, 15.0]
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
    while length(p_bootstrap) < n
        Threads.@threads for i=1:n-length(p_bootstrap)
            boot_df = block_bootstrap(df, 20, 5)
            try
                myfit = fit_data(boot_df.time, boot_df.avg_height)
                append!(p_bootstrap, [myfit.u])
                print(length(p_bootstrap))
            catch
            end
        end
    end
    return reduce(hcat, p_bootstrap)'
end
strain_name = "bh1514"
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== strain_name && x.time .< 48 &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
#@df df scatter(:time, :avg_height, group=(:replicate), legend=false, xlabel="Time [hr]",
#               ylabel="Height [μm]", size=(500, 400))
data = boot_fit(df, 10)

##
parameters_frame = DataFrame("α"=>data[:,1], "β"=>data[:,2], "L"=>data[:,3])
parameters_frame.h_max = parameters_frame.α .* parameters_frame.L ./ parameters_frame.β
CSV.write("data/sims/bootstrap/boot_"*strain_name*".csv", parameters_frame)

#= 
pf = DataFrame("strain"=>String[],"name"=>String[], "α"=>Float64[],
               "β"=>Float64[], "L"=>Float64[], "h_max"=>Float64[]) 
push!(pf, ["bgt127" "fit_48" bgt127_48[1] bgt127_48[2] bgt127_48[3] hmax(bgt127_48)])
push!(pf, ["bgt127" "fit_all" bgt127_all[1] bgt127_all[2] bgt127_all[3] hmax(bgt127_all)])
push!(pf, ["jt305" "fit_48" jt305_48[1] jt305_48[2] jt305_48[3] hmax(jt305_48)])
push!(pf, ["jt305" "fit_all" jt305_all[1] jt305_all[2] jt305_all[3] hmax(jt305_all)])
push!(pf, ["gob33" "fit_48" gob33_48[1] gob33_48[2] gob33_48[3] hmax(gob33_48)])
push!(pf, ["gob33" "fit_all" gob33_all[1] gob33_all[2] gob33_all[3] hmax(gob33_all)])
data_bgt127 = DataFrame(CSV.File("data/sims/bootstrap/boot_bgt127.csv"))
data_bgt127.strain = repeat(["bgt127"], 1000)
numbers = [lpad(i, 4, "0") for i=1:1000]
data_bgt127.name = numbers
append!(pf, data_bgt127)
data_jt305 = DataFrame(CSV.File("data/sims/bootstrap/boot_jt305.csv"))
data_jt305.strain = repeat(["jt305"], 1000)
data_jt305.name = numbers
append!(pf, data_jt305)
data_gob33 = DataFrame(CSV.File("data/sims/bootstrap/boot_gob33.csv"))
data_gob33.strain = repeat(["gob33"], 1000)
data_gob33.name = numbers
append!(pf, data_gob33)

CSV.write("data/sims/bootstrap/allfits.csv", pf) =#
##

