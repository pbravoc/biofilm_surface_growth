using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Plots, StatsPlots 

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
function interface_n(du, u, p, t)
    h, c = u 
    α, β, hstar, K_c, ϵ = p
    du[1] = α*G.(h, hstar)*(c/(K_c + c)) - β*h 
    du[2] = -ϵ*(c/(K_c + c)) 
    return du
end
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function fit_data(t_data, h_data, model, pguess=[0.8, 0.1, 15.0])
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
    result_ode = DiffEqFlux.sciml_train(loss, pguess,
                        lower_bounds = [0.1, 50.0],
                        upper_bounds = [5.0, 500.0])
    return result_ode
end 
##
df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x-> x.strain .== "bgt127" && x.time .< 48 &&
                x.replicate in ["A", "B", "C"], df)
##
plot()
for i=1:20
    @df df[[i, 100+i, 199+i], :] scatter!(:time, :avg_height, color=:black, label=false, alpha=0.2)
    plot!(xlim=(0.0, 48.0), ylim=(-1, 220.0), xlabel="Time [hr]", ylabel="Height [μm]")
end
anim = @animate for i=21:98
    @df df[[i, 100+i, 199+i], :] scatter!(:time, :avg_height, color=:black, label=false, alpha=0.2)
end
gif(anim, "figs/animations/growth_gif.gif", fps = 10)
##
df_reference = filter(x-> x.replicate in ["C"], df)
my_pars2 = []
for i=25:98
    tf = filter(x-> x.time .<= df_reference.time[i], df)
    fit_pars = fit_data(tf.time, tf.avg_height, logistic, [0.1, 200.0])
    #fit_pars = fit_data(tf.time, tf.avg_height, interface)
    print(fit_pars)
    append!(my_pars2, [fit_pars])
end
##
s = 5
anim = @animate for i=s:length(my_pars)
    t_end = df_reference.time[25+i]
    tf = filter(x-> x.time .<= t_end, df)
    @df tf scatter(:time, :avg_height, color=:black, label=false, alpha=0.2, xlim=(-0.5, 48.0), ylim=(0.0, 220.0))
    prob = ODEProblem(interface, [0.1], (0.0, 48.0), my_pars[i])
    sol = solve(prob, saveat=0.5, save_idxs=1)
    idxs = sol.t .<= t_end
    idxs2 = sol.t .> t_end
    #plot!(sol.t[idxs], sol.u[idxs], color=:red, label=false, title=i,linewidth=2)
    #plot!(sol.t[idxs2], sol.u[idxs2], color=:red, label=false, linestyle=:dash, linewidth=2)
    prob = ODEProblem(logistic, [0.1], (0.0, 48.0), my_pars2[i])
    sol = solve(prob, saveat=0.5, save_idxs=1)
    plot!(sol.t[idxs], sol.u[idxs], color=:blue, label=false, title=i, linewidth=2)
    plot!(sol.t[idxs2], sol.u[idxs2], color=:blue, label=false, linestyle=:dash, linewidth=2)
    plot!(xlabel="Time [hr]", ylabel="Height[μm]", title=string(round(t_end, digits=1))*" hours", size=(500, 400), dpi=500)
    end
gif(anim, "figs/animations/growth_logistic_fit.gif", fps = 10)
##
tf = filter(x-> x.time .<= df_reference.time[30], df)
@df tf scatter(:time, :avg_height, color=:black, label=false, alpha=0.2)
prob = ODEProblem(logistic, [0.1], (0.0, 48.0), my_pars2[5])
sol = solve(prob, saveat=0.1)
#plot!(sol, color=:blue, label=false, title=30, linewidth=2)
prob = ODEProblem(interface, [0.1], (0.0, 11.0), my_pars[5])
sol = solve(prob, saveat=0.1)
#plot!(sol, color=:red, label=false, title=i, linewidth=2)
plot!(xlabel="Time [hr]", ylabel="Height[μm]", title="Early times", size=(500, 400), xlim=(-0.5, 10.5), ylim=(-0.5, 58.0))
#savefig("figs/animations/pre_growth0.svg")