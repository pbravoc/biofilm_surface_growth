using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using Plots, StatsPlots
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
    result_ode = [NaN, NaN, NaN]
    u0 = [0.2]
    pguess = [0.8, 0.1, 15.0]
    prob = ODEProblem(interface, u0, (0.0, t_data[end]), pguess)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, pguess, 
                                        lower_bounds = [1e-4, 1e-4, 1e-4], 
                                        upper_bounds = [10.0, 2.0, 1000.0])
    return result_ode
end 
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== "jt305" && x.time .< 48 &&
                x.replicate in ["A", "B", "C"], Df)
##
@df df scatter(:time, :avg_height, group=(:replicate), legend=false, xlabel="Time [hr]",
               ylabel="Height [μm]", size=(500, 400), xlim=(-1, 49), ylim=(-1, 125))
savefig("figs/bootstrapping/data_sample.svg")
##
anim = @animate for i=1:30
    b_data = block_bootstrap(df, 20, 5)
    @df b_data scatter(:time, :avg_height, group=(:replicate), legend=false, xlabel="Time [hr]",
                        ylabel="Height [μm]", size=(500, 400), xlim=(-1, 49), ylim=(-1, 125),dpi=300)
end
gif(anim, "figs/bootstrapping/sample_data.gif", fps = 5)
##
params = fit_data(b_data.time, b_data.avg_height)
##
params = fit_data(df.time, df.avg_height)
##
N = 1000
Parameters = []
##
for i = 1:N
    b_data = block_bootstrap(df, 20, 5)
    #@df b_data scatter(:time, :avg_height, group=(:replicate), legend=false)
    curr_params = fit_data(b_data.time, b_data.avg_height)
    append!(Parameters, [curr_params.u])
    print(i)
end
##
p1 = histogram(data[1,:], title="α")
p2 = histogram(data[2,:], title="β")
p3 = histogram(data[3,:], title="L") 
p4 = histogram(parameters_frame.h_max, title="h_max") 
plot(p1, p2, p3, p4)
#savefig("figs/bayesian_fitting/boot_gob33.svg")

##
data = reduce(hcat, Parameters)
parameters_frame = DataFrame("α"=>data[1,1:1000], "β"=>data[2,1:1000], "L"=>data[3,1:1000])
parameters_frame.h_max = parameters_frame.α .* parameters_frame.L ./ parameters_frame.β
CSV.write("data/sims/bootstrap/boot_gob33.csv", parameters_frame)
##
@df parameters_frame scatter(:L, :h_max, alpha=0.3)
##
@df parameters_frame marginalkde(:α, :β)
##
p4 = histogram(parameters_frame.h_max, xlim=(0, 2000)) 
#savefig("figs/bayesian_fitting/bootzoom_gob33.svg")

##
#gob33_48 = params
#bgt127_48 = params
#jt305_48 = params
#gob33_all = params
#bgt127_all = params
jt305_all = params
##
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
##
CSV.write("data/sims/bootstrap/allfits.csv", pf)
##

