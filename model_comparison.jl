#=
This script generates a new dataset that contains 
the avg + std height from the aggregated data (A, B, C), and 
the respective best-fit predictions for multiple models in the 48h
range.
=# 

using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors
using LsqFit
using DifferentialEquations

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
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
model_list = [interface]
df = DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.time .< 48 && x.avg_height .> 0 && 
                x.replicate in ["A", "B", "C"], df)  
pf = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
Data = DataFrame(strain = String[], time=Float32[], time_error=Float32[],
                 avg_height = Float32[], std_height = Float32[], 
                 interface = Float32[]) # Add the other models
##
for strain in unique(df.strain)
    print(strain)
    tf = filter(x->x.strain .== strain, df)
    p = Array(filter(x->x.fit.=="48h" && x.strain .== strain, 
                     pf)[:,3:5])
    t, t_e, h, h_e = order_average(tf)
    sf = DataFrame()
    sf.strain = repeat([strain], length(t))
    sf.time = t 
    sf.time_error = t_e 
    sf.avg_height = h 
    sf.std_height = h_e
    problem = ODEProblem(interface, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.interface = sol.u
    global Data = vcat(Data, sf)
end
##
CSV.write("data/timelapses/model_predictions.csv", Data)
