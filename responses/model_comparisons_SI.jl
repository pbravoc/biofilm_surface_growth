using DataFrames, CSV
using Statistics, NaNMath
using DifferentialEquations

G(z, zstar) = z < zstar ? z : zstar 
function nutrient_n(du, u , p, t)
    h, c = u
    α, β, K_c, ϵ = p         
    du[1] = α*h*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end
function logistic_n(du, u , p, t)
    h, c = u
    α, K_h, K_c, ϵ = p         
    du[1] = α*h*(1- h/(K_h))*(c/(K_c + c)) 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du 
end
function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    #du[1] = α * h - (α*h*h/K_h + h)
    du[1] = α*h*(1- h/(K_h))
    #du[1] = α * h *(1- (h/(K_h + h)))
    return du 
end

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
df = DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.time .< 48 && x.avg_height .> 0 && 
                x.replicate in ["A", "B", "C"], df)  
params_interface = DataFrame(CSV.File("data/timelapses/parameters/params_interface.csv"))
params_logisticnd = DataFrame(CSV.File("data/timelapses/parameters/params_logisticnd.csv"))
params_logistic = DataFrame(CSV.File("data/timelapses/parameters/params_logistic.csv"))
params_nutrient = DataFrame(CSV.File("data/timelapses/parameters/params_nutrient.csv"))

Data = DataFrame(strain = String[], time=Float32[], time_error=Float32[],
                 avg_height = Float32[], std_height = Float32[], 
                 interface = Float32[], logistic_nd = Float32[],
                 logistic = Float32[], nutrient = Float32[]) 
##
for strain in unique(df.strain)
    sf = DataFrame()
    println(strain)
    # Get average values
    tf = filter(x->x.strain .== strain, df)
    t, t_e, h, h_e = order_average(tf)
    sf.strain = repeat([strain], length(t))
    sf.time = t 
    sf.time_error = t_e 
    sf.avg_height = h 
    sf.std_height = h_e

    # Interface
    p = Array(filter(x->x.strain .== strain, params_interface)[:,3:end])
    problem = ODEProblem(interface, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.interface = sol.u

    # Logistic + Nutrient depletion
    p = Array(filter(x->x.strain .== strain, params_logisticnd)[:,3:end])
    problem = ODEProblem(logistic_n, [0.1, 1.0], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.logistic_nd = sol.u

    # Logistic
    p = Array(filter(x->x.strain .== strain, params_logistic)[:,3:end])
    problem = ODEProblem(logistic, [0.1], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.logistic = sol.u

    # Nutrient depletion
    p = Array(filter(x->x.strain .== strain, params_nutrient)[:,3:end])
    problem = ODEProblem(nutrient_n, [0.1, 1.0], (0.0, t[end]), p)
    sol = solve(problem, saveat=t, save_idxs=1)
    sf.nutrient = sol.u

    append!(Data, sf)
end
##
CSV.write("data/timelapses/model_predictions_all.csv", Data)
