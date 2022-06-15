#=
This script generates a new dataset that contains 
the avg + std height from the aggregated data (A, B, C), and 
the respective best-fit predictions for multiple models in the 48h
range.

For bgt127, jt305, and gob33 we also calculate using the ALL fit.
=# 

using DataFrames, CSV
using Statistics, NaNMath
using DifferentialEquations

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    du[1] = α*h*(1- h/(K_h))
    return du 
end

"Gets the average depending on the order, the advantage of this 
is that if one replicate is void, then that point will not be 
taken in account for the calculation. "
function order_average(Df)
    n = maximum(Df.order)
    t, t_e, h, h_e = zeros(n), zeros(n), zeros(n), zeros(n)
    for i=1:maximum(Df.order)
        tf = filter(x->x.order .== i, Df)
        t[i], t_e[i] = mean(tf.time), std(tf.time)
        h[i], h_e[i] = mean(tf.avg_height), std(tf.avg_height)
    end
    return t, t_e, h, h_e
end
##
longtime_list = ["bgt127", "gob33", "jt305"]
df = DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.time .< 48 && x.avg_height .> 0 && 
                x.replicate in ["A", "B", "C"], df)  
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
pf2 = DataFrame(CSV.File("data/timelapses/fit_params_logistic.csv"))
Data = DataFrame(strain = String[], time=Float32[], time_error=Float32[],
                 avg_height = Float32[], std_height = Float32[], 
                 interface = Float32[], interface_long = Float32[],
                 interface_alt = Float32[],
                 logistic = Float32[], logistic_long = Float32[]) 
##
for strain in unique(df.strain)
    sf = DataFrame()
    print(strain)
    # Get average values
    tf = filter(x->x.strain .== strain, df)
    t, t_e, h, h_e = order_average(tf)
    sf.strain = repeat([strain], length(t))
    sf.time = t 
    sf.time_error = t_e 
    sf.avg_height = h 
    sf.std_height = h_e
    sf.interface_long = zeros(length(t)) .+= NaN
    sf.logistic_long = zeros(length(t)) .+= NaN
    # Get 48h fits
    p = Array(filter(x->x.fit.=="48h" && x.strain .== strain, 
                     pf)[:,3:5])
    problem = ODEProblem(interface, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.interface = sol.u
    # Interface alt 
    p = Array(filter(x->x.fit.=="alt" && x.strain .== strain, 
                     pf)[:,3:5])
    problem = ODEProblem(interface, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.interface_alt = sol.u
    ## Logistic
    p = Array(filter(x->x.fit.=="48h" && x.strain .== strain, 
                     pf2)[:,3:4])
    problem = ODEProblem(logistic, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.logistic = sol.u

    ## Get long fits
    if strain in longtime_list
        p = Array(filter(x->x.fit.=="long" && x.strain .== strain, 
                        pf)[:,3:5])
        problem = ODEProblem(interface, [0.1], (0.0, t[end]), p)
        sol = solve(problem, saveat=t, save_idxs=1)
        sf.interface_long = sol.u
        p = Array(filter(x->x.fit.=="long" && x.strain .== strain, 
                         pf2)[:,3:4])
        problem = ODEProblem(logistic, [0.1], (0.0, t[end]), p)
        sol = solve(problem, saveat=t, save_idxs=1)
        sf.logistic_long = sol.u
    end
    append!(Data, sf)
end
##
CSV.write("data/timelapses/model_predictionsv2.csv", Data)
