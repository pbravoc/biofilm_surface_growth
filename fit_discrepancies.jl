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
function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    du[1] = α*h*(1- h/(K_h))
    return du 
end
function get_prediction(model, p, t)
    prob = ODEProblem(model, [0.1], (0.0, 48.0), p)
    sol = solve(prob, saveat=t, save_idxs=1)
    return sol.u
end
function order_average(Df, strain_name)
    df = filter(x->x.strain .== strain_name,Df)
    n = maximum(df.order)
    t, t_e, h, h_e = zeros(n), zeros(n), zeros(n), zeros(n)
    for i=1:maximum(df.order)
        tf = filter(x->x.order .== i,df)
        t[i], t_e[i] = mean(tf.time), std(tf.time)
        h[i], h_e[i] = mean(tf.avg_height), std(tf.avg_height)
    end
    return t, h, h_e
end

df = DataFrame(CSV.File("data/timelapses/database.csv"))
pf_logistic = DataFrame(CSV.File("data/timelapses/fit_params_logistic.csv"))
pf_interface = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
order = []
for strain in unique(df.strain)
    for repli in unique(df.replicate)
        tf = filter(row->row.replicate.==repli && row.strain .== strain, df);
        my_order = Array(1:size(tf)[1])
        append!(order, my_order)  
    end  
end
df.order = order
df = filter(x-> x.time .< 48 && x.avg_height .> 0 && 
                x.replicate in ["A", "B", "C"], df)  
summary_df = DataFrame("strain"=>String[], "t"=>Float64[], 
                       "h"=>Float64[], "h_std"=>Float64[], 
                       "interface_all"=>Float64[],
                       "interface_48"=>Float64[],
                       "logistic_all"=>Float64[],
                       "logistic_48"=>Float64[])
##                       
for strain in unique(df.strain)
    tf = filter(x-> x.strain .== strain, df) 
    t, h, h_e = order_average(tf, strain)
    data_list = []
    for fit_type in ["all", "48h"]
        p_in = filter(x->x.strain .== strain && x.fit .== fit_type, 
                    pf_interface)[1,3:5]
        p_log = filter(x->x.strain .== strain && x.fit .== fit_type, 
                    pf_logistic)[1,3:4]
        append!(data_list, [get_prediction(interface, p_in, t)])
        append!(data_list, [get_prediction(logistic, p_log, t)])
    end 
    tf = DataFrame("strain"=>repeat([strain], length(t)), "t"=>t, 
                    "h"=>h, "h_std"=>h_e, 
                    "interface_all"=>data_list[1],
                    "interface_48"=>data_list[3],
                    "logistic_all"=>data_list[2],
                    "logistic_48"=>data_list[4])
    append!(summary_df, tf)
end

##
#CSV.write("data/timelapses/fit_predictions.csv", summary_df)
## Long time fits
function get_prediction_long(model, p)
    prob = ODEProblem(model, [0.1], (0.0, 350.0), p)
    sol = solve(prob, saveat=0.5, save_idxs=1)
    return sol.u
end
summary_df = DataFrame("strain"=>String[], "t"=>Float64[], 
                       "interface_all"=>Float64[],
                       "interface_48"=>Float64[],
                       "logistic_all"=>Float64[],
                       "logistic_48"=>Float64[])
##                       
for strain in unique(pf_interface.strain)
    data_list = []
    for fit_type in ["all", "48h"]
        p_in = filter(x->x.strain .== strain && x.fit .== fit_type, 
                    pf_interface)[1,3:5]
        p_log = filter(x->x.strain .== strain && x.fit .== fit_type, 
                    pf_logistic)[1,3:4]
        append!(data_list, [get_prediction_long(interface, p_in)])
        append!(data_list, [get_prediction_long(logistic, p_log)])
    end 
    tf = DataFrame("strain"=>repeat([strain], length(Array(0:0.5:350.0))), "t"=>Array(0:0.5:350.0), 
                    "interface_all"=>data_list[1],
                    "interface_48"=>data_list[3],
                    "logistic_all"=>data_list[2],
                    "logistic_48"=>data_list[4])
    append!(summary_df, tf)
end
##
CSV.write("data/timelapses/longfit_predictions.csv", summary_df)
