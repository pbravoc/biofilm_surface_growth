#=
We calculate the trajectories for the bootstrapped data for 
bgt127, jt305, and gob33. Quantile calculation is best left 
for the plotting part of the pipeline. 
48h best fit, and ALL fit also calculated here.
=#
using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob

G(z, zstar) = z < zstar ? z : zstar 
function interface_boot(du, u, p, t)
    h = u 
    α, β, hstar = p[:,1], p[:,2], p[:,3]
    for i=1:1000 # I know this is bad... but it doesnt work otherwise lol!
        du[i] = α[i] .* G.(h[i], hstar[i]) .- β[i] .* h[i] 
    end
    return du
end

function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
##
pf2 = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
df_long = DataFrame(CSV.File("data/timelapses/longtime_data.csv")) # Experimental
strain_list = ["bgt127", "jt305", "gob33"]
Data = DataFrame()
u0 = zeros(1000) .+ 0.1
for strain in strain_list
    # bootstrapped trajectories 
    t_save = df_long[df_long.strain .== strain,:].time
    #t_save = 0.5
    pf = DataFrame(CSV.File("data/sims/bootstrap/boot_"*strain*".csv"))
    tf = filter(x->x.strain .== strain, pf2)
    p = Array(pf)
    problem = ODEProblem(interface_boot, u0, (0.0, 350.0), p)
    sol = solve(problem, saveat=t_save)
    t = sol.t
    sol = reduce(hcat, sol.u)
    number_list = "id_" .*lpad.(string.(Array(0:1:999)), 3, '0')
    df= DataFrame(Dict([number_list[i]=>sol[i,:] for i=1:1000]))
    # 48h best fit trajectories 
    p = Array(tf[tf.fit .== "48h", 3:5])
    problem = ODEProblem(interface, u0, (0.0, 350.0), p)
    sol = solve(problem, saveat=t_save, save_idxs=1)
    df.t48 = sol.u
    # all best fit trajectories 
    p = Array(tf[tf.fit .== "long", 3:5])
    problem = ODEProblem(interface, u0, (0.0, 350.0), p)
    sol = solve(problem, saveat=t_save, save_idxs=1)
    df.tall = sol.u
    df.time = t 
    df.strain = repeat([strain], length(t))
    append!(Data, df)
end
##
#CSV.write("data/sims/bootstrap/boot_trajectories.csv", Data)
CSV.write("data/sims/bootstrap/boot_trajectories_long_ref.csv", Data)

##
using Plots, StatsPlots 
@df Data plot(:time, [:t48, :tall], group=:strain, color=:auto)
##
