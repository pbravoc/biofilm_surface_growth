using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes
using DifferentialEquations

function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    return t, h_avg, h_std
end
function logistic(du, u , p, t)
    h = u[1]
    α, K_h = p  
    #du[1] = α * h - (α*h*h/K_h + h)
    du[1] = α*h*(1- h/(K_h))
    #du[1] = α * h *(1- (h/(K_h + h)))
    return du 
end
G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.replicate in ["A", "B", "C"] && 
               x.time .<= 48, Df)

##
data_bgt127 = get_average(df, "bgt127")
data_gob33 = get_average(df, "gob33")
data_jt305 = get_average(df, "jt305")

plot()
for data in [data_bgt127, data_gob33, data_jt305]
    scatter!(data[1], data[2], yerror=data[3])
end
plot!(xlabel="Time [hr]", ylabel="Height [μm]")

## 48h fit
parse_params(file, i, j) = parse.(Float64, split(file[i,j][2:end-1], ", "))
# bgt127
function fit_48(df, strain_name)
    data = get_average(df, strain_name)
    params = DataFrame(CSV.File("data/sims/"*strain_name*"_params.csv"))
    log_pars = parse_params(params, 4,2) 
    in_pars = parse_params(params, 5,2) 
    u0 = [0.2]
    log_problem = ODEProblem(logistic, u0, [0.0, 350.0], log_pars)
    in_problem = ODEProblem(interface, u0, [0.0, 350.0], in_pars)
    sol_log = solve(log_problem, saveat=data[1], save_idxs=1)
    sol_log_long = solve(log_problem, saveat=1.0, save_idxs=1)
    sol_in = solve(in_problem, saveat=data[1], save_idxs=1)
    sol_in_log = solve(in_problem, saveat=1.0, save_idxs=1)
    return data[1], data[2], data[3], sol_log.u, sol_in.u, 
           sol_log_long, sol_in_log
end
function plot_48(data, strain_name)
    tf = filter(x->x.strain .== strain_name, Df)
    p1 = @df tf scatter(:time, :avg_height, legend=:bottomright, markersize=3,
                   color=:black, alpha=0.5, label=false)
    plot!(data[7], color=1, linewidth=2, label="Interface")
    plot!(data[6], color=2, linewidth=2, label="Logistic")
    plot!(legend=:bottomright, title=strain_name, xlabel="Time [hr]", ylabel= "Height [μm]")
    p2 = scatter(data[1], data[2]-data[2], yerror=data[3], markersize=0, label=false)
    plot!(data[1], data[5]-data[2], color=1, linewidth=2, xlabel="Time [hr]")
    plot!(data[1], data[4]-data[2], color=2, linewidth=2, ylabel="Residual [μm]",  legend=false)
    return plot(p1, p2, layout=(2,1), size=(400, 500))
end

bgt127_data = fit_48(df, "bgt127")
jt305_data = fit_48(df, "jt305")
gob33_data = fit_48(df, "gob33")

plot(plot_48(bgt127_data, "bgt127"), 
     plot_48(gob33_data, "gob33"), 
     plot_48(jt305_data, "jt305"), 
     layout=(1,3), size=(850, 500))
#savefig("figs/longfit_48.svg")
##


##
function fit_all(df, strain_name)
    data = get_average(df, strain_name)
    log_all = DataFrame(CSV.File("data/sims/allpointfit_log.csv"))
    in_all = DataFrame(CSV.File("data/sims/allpointfit.csv"))
    log_pars = Array(log_all[log_all.Name .== strain_name,2:3 ])
    in_pars = Array(in_all[in_all.Name .== strain_name,2:4 ])
    u0 = [0.2]
    log_problem = ODEProblem(logistic, u0, [0.0, 350.0], log_pars)
    in_problem = ODEProblem(interface, u0, [0.0, 350.0], in_pars)
    sol_log = solve(log_problem, saveat=data[1], save_idxs=1)
    sol_log_long = solve(log_problem, saveat=1.0, save_idxs=1)
    sol_in = solve(in_problem, saveat=data[1], save_idxs=1)
    sol_in_log = solve(in_problem, saveat=1.0, save_idxs=1)
    return data[1], data[2], data[3], sol_log.u, sol_in.u, 
           sol_log_long, sol_in_log
end

bgt127_data = fit_all(df, "bgt127")
gob33_data = fit_all(df, "gob33")
jt305_data = fit_all(df, "jt305")

plot(plot_48(bgt127_data, "bgt127"), 
     plot_48(gob33_data, "gob33"), 
     plot_48(jt305_data, "jt305"), 
     layout=(1,3), size=(850, 500))
savefig("figs/longfit_all.svg")
##

function residual_compare(strain_name)
    data = fit_48(df, strain_name)
    data_all = fit_all(df, strain_name)
    p1 = scatter(data[1], data[2]-data[2], yerror=data[3], markersize=0, label=false)
    plot!(data[1], data[5]-data[2], color=1, linewidth=2, xlabel="Time [hr]", label="Interface_48",title=strain_name)
    plot!(data[1], data[4]-data[2], color=2, linewidth=2, ylabel="Residual [μm]", label="Logistic_48",legend=:topleft)
    plot!(data_all[1], data_all[5]-data_all[2], color=1, linestyle=:dash, linewidth=2,label="Interface",xlabel="Time [hr]")
    plot!(data_all[1], data_all[4]-data_all[2], color=2, linestyle=:dash, linewidth=2,label="Logistic", xlabel="Time [hr]", ylim=(-20,20))
    return p1
end

plot(residual_compare("bgt127"), 
     residual_compare("gob33"), 
     residual_compare("jt305"), 
     layout=(1,3), size=(1000, 400), left_margin=4mm, bottom_margin=5mm)
savefig("figs/residuals_constrained.svg")
