using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots
using DifferentialEquations, DiffEqFlux

function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    return t, h_avg, h_std
end
G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function fit_data(t_data, h_data)
    result_ode = [NaN, NaN, NaN]
    try
        u0 = [0.2]
        pguess = [0.8, 0.07, 15]
        prob = ODEProblem(interface, u0, (0.0, t_data[end]), pguess)
        function loss(p)
            sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
            sol_array = reduce(vcat, sol.u)
            loss = sum(abs2, sol_array .- h_data)
            return loss, sol
        end
        result_ode = DiffEqFlux.sciml_train(loss, pguess,  
                                            lower_bounds = [0.0, 0.0, 0.0], 
                                            upper_bounds = [1.0, 0.5, 50.0])
    catch e 
        print("poto")
    end
    return result_ode
end 
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.replicate in ["A", "B", "C"] && x.strain .== "jt305", Df)
t, h, dh = get_average(df, "jt305")
param_predictions = zeros(length(t), 3)
scatter(t, h, yerror=dh)

##
2+2
##
Threads.@threads for i = 10:length(t)
    print(i)
    my_params = fit_data(t[1:i], h[1:i])
    param_predictions[i,:] = my_params
    print(my_params)
end
##
h_eq = param_predictions[:,1] .* param_predictions[:,3] ./ param_predictions[:,2]
scatter(h_eq[10:50])
##
pf = DataFrame("time"=>t, "data"=>h, "α"=>param_predictions[:,1],
               "β"=>param_predictions[:,2], "L"=>param_predictions[:,3])
CSV.write("data/sims/tp_bounded.csv", pf)


##
pf 