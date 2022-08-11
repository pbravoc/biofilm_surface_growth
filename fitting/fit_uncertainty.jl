using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob

function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size-1) for x in sample_start])
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
function fit_data(t_data, h_data, pguess=[0.8447392179913538,0.062439602238195,15.306008084305928])
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
    result_ode = DiffEqFlux.sciml_train(loss, pguess,
                                        maxiters=100)
    return result_ode
end 

function boot_fit(df, n, n_blocks)
    p_bootstrap = []
    while length(p_bootstrap) < n
        Threads.@threads for i=1:n-length(p_bootstrap)
            boot_df = block_bootstrap(df, n_blocks, 5)
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
## This is to do bootstrap fitting on 48h data
strain_name = "bgt127"
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== strain_name &&
                x.replicate in ["A", "B", "C"] &&
                x.avg_height > 0 , Df)
##
n_boots = 100
P = []
PF = DataFrame()
for t in 10:2:48
    tf = filter(x-> x.time < t , df)
    n_target = Int(floor(size(tf)[1]/3))
    n_blocks = Int(floor(n_target/5))
    p_bootstrap = []      
    while length(p_bootstrap) < n_boots
        Threads.@threads for i=1:n_boots-length(p_bootstrap)
            try
                boot_df = block_bootstrap(tf, n_blocks, 5)
                myfit = fit_data(boot_df.time, boot_df.avg_height)
                append!(p_bootstrap, [myfit.u])
            catch
            end
        end
    end 
    p_bootstrap = reduce(hcat, p_bootstrap)'
    pf = DataFrame("t"=>repeat([t], n_boots), 
                   "h"=>repeat([maximum(tf.avg_height)], n_boots), 
                   "α"=>p_bootstrap[:,1], 
                   "β"=>p_bootstrap[:,2], 
                   "L"=>p_bootstrap[:,3])
    global PF = append!(PF, pf)
end
CSV.write("data/sims/bootstrap/timeboot_bgt127.csv", PF)
##
using Plots, StatsPlots 
using Statistics

function summary_long(df, strain_name, L, target)
    df.h_max = df.α .* df.L ./ df.β
    h_low = []
    h_mean = []
    h_high = []
    t_sample = []
    h_sample = []
    h_sample_L = []
    h_sample_P = []
    pred_error = []
    pred_error_P = []
    h_error = []
    h_error_P = []
    for t in unique(df.t)
        tf = filter(x-> x.t .== t, df)
        h_array = tf.h_max
        h = sort(h_array)[30:70]
        hl, hh = quantile(h_array, [0.3, 0.7])
        append!(t_sample, tf.t[1])
        append!(h_sample, tf.h[1])
        append!(h_mean, mean(h))
        append!(h_low, hl)
        append!(h_high, hh)
        append!(pred_error, mean(h)-target)
        append!(pred_error_P, abs.(mean(h) - target)/target)
        append!(h_error, std(h))
        append!(h_error_P, std(h)/target)
    end
    pf = DataFrame("strain"=>repeat([strain_name], length(t_sample)),
                   "t_sample"=>t_sample, "h_sample"=>h_sample, 
                   "h_low"=>h_low, "h_mean"=>h_mean, 
                   "h_high"=>h_high,
                   "err_low"=>abs.(h_mean-h_low), 
                   "err_high"=>abs.(h_mean-h_high),
                   "h_sample_L"=>h_sample ./ L, "h_sample_P"=>h_sample ./ target,
                   "pred_error"=>pred_error, "pred_error_P"=>pred_error_P,
                   "h_error"=>h_error, "h_error_P"=>h_error_P)
    return pf 
end
##
strain_list = ["bgt127", "jt305", "gob33"]
L_list = [15.30, 12.13, 28.74]
target_list = [207.7, 322.8, 831.8]
data = DataFrame()
for i=1:3
    temp_frame = DataFrame(CSV.File("data/sims/bootstrap/timeboot_"*strain_list[i]*".csv"))
    pf = summary_long(temp_frame, strain_list[i], L_list[i], target_list[i])
    append!(data, pf)
end
##
mc = [ColorSchemes.okabe_ito[5], ColorSchemes.okabe_ito[6], ColorSchemes.okabe_ito[4]]'
ms = [:diamond :square :circle]
ml = ["Aeromonas" "Yeast (aa)" "E coli"]
plot()
p1 = @df data scatter(:h_sample_L, :pred_error_P, group=:strain,
                 color=mc, marker=ms, label=ml, xlabel="L", ylabel="PT error (%)")
p2 = @df data scatter(:h_sample_L, :h_error_P, group=:strain,
                 color=mc, marker=ms, label=ml, xlabel="L", ylabel="P range (%)")
p3 = @df data scatter(:h_sample_P, :pred_error_P, group=:strain,
                 color=mc, marker=ms, label=ml, xlabel="h_max", ylabel="PT error (%)")
p4 = @df data scatter(:h_sample_P, :h_error_P, group=:strain,
                 color=mc, marker=ms, label=ml, xlabel="h_max", ylabel="P range (%)")
plot(p1, p2, p3, p4, legend=false, dpi=500)
savefig("figs/fig4/pt_all.png")
##
@df data scatter(:pred_error_P, :h_error_P, group=:strain,
                 color=mc, marker=ms, label=ml, ylabel="P range", xlabel="PT error (%)",
                 xscale=:log, yscale=:log)