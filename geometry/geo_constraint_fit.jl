using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob
using Plots, SpecialFunctions

function c(x, t)
    L = sqrt(D*t)
    return erfc.(x / L)
end

function int_c(h, L)
    t1 = L*(-exp.(-h^2 / (L^2)))
    t2 = sqrt(π) * h * erfc(h/L)
    return (t1 + t2 + L)
end

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

function interface_full(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*int_c.(h, hstar) - β*h 
    return du
end

function fit_data(t_data, h_data, model, pguess=[0.8, 0.05, 15.0])
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

##
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0, Df)     
df = filter(x-> x.strain .== "bgt127" && x.time .<48, Df)

##
fit_interface = fit_data(df.time, df.avg_height, interface)
##
fit_interfacef = fit_data(df.time, df.avg_height, interface_full)
##
prob1 = ODEProblem(interface, u0, (0.0, 48.0), fit_interface)
sol1 = solve(prob1, save_idxs=1, saveat=0.05)
prob2 = ODEProblem(interface_full, u0, (0.0, 48.0), fit_interfacef)
sol2 = solve(prob2, save_idxs=1, saveat=0.05)
##
p1 = plot(sol2, linewidth=2.0, alpha=1.0, legend=:topleft, color=:black, label="Analyitical")
plot!(sol1, linewidth=2.0, alpha=1.0, color=:red, linestyle=:dash, label="Approximation")
plot!(xlabel="Time [hr]", ylabel="Height [μm]")
lens!([1, 7], [0, 20], inset = (1, bbox(0.55, 0.4, 0.4, 0.5)), subplot = 2)
p2 = plot(sol2.t, sol2.u-sol1.u, color=:black, linewidth=2.0, label=false)
hline!([0.0], color=:gray, linestyle=:dash, linewidth=1.5, label=false)
plot!(xlabel="Time [hr]", ylabel="Difference [μm]")
plot(p1, p2, grid=false, size=(700, 300), bottom_margin=3mm, left_margin=4mm)
savefig("figs/fig3/geo_constraint.pdf")