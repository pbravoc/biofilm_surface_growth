using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots
using DifferentialEquations, DiffEqFlux

G(z, zstar) = z < zstar ? z : zstar 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end
function get_average(df, strain_name, n)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/n)
    h = reshape(abs.(tf.avg_height), (l, n))
    print(minimum(h))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    return t, h_avg, h_std
end
function fit_data(t_data, h_data)
    result_ode = [NaN, NaN, NaN]
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], abs.(h_data[idxs])  # This is needed to compare 1-1
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
    u0 = [0.2]
    prob = ODEProblem(interface, u0, (0.0, t_data[end]), result_ode)
    sol = solve(prob, Tsit5(), saveat=0.5, save_idxs=1)
    return result_ode, sol
end 

Df =  filter(x->x.time .<= 47.97 && x.replicate in ["A", "B", "C"], DataFrame(CSV.File("data/timelapses/database.csv")))
Data = []
sols = []
##
df = filter(x->x.strain .== "bacillus" && 
               x.replicate in ["B", "C"], Df)
t, h, h_e = get_average(df, "bacillus", 2)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])

df = filter(x->x.strain .== "gob33" , Df)
t, h, h_e = get_average(df, "gob33", 3)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])

df = filter(x->x.strain .== "jt1080" , Df)
t, h, h_e = get_average(df, "jt1080", 3)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])

df = filter(x->x.strain .== "sn503"  && 
            x.replicate in ["A", "B"], Df)
df[109, :avg_height] = 33.0
t, h, h_e = get_average(df, "sn503", 2)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])

df = filter(x->x.strain .== "jt305" , Df)
t, h, h_e = get_average(df, "jt305", 3)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])

df = filter(x->x.strain .== "bgt127" , Df)
t, h, h_e = get_average(df, "bgt127", 3)
p, s = fit_data(df.time, df.avg_height)
append!(Data, [[t, h, h_e]])
append!(sols, [s])
##
plot()
plot!(sols[1])
plot!(sols[2])
plot!(sols[3])
plot!(sols[4])
plot!(sols[5])
plot!(sols[6])
plot!(xlabel="Time [hr]", ylabel="Height [μm]")
##
my_plots = []
strain_names = ["Bacillus", "Yeast (aa)", "Vcholerae (EPS-)", "Vcholerae (wt)", "Ecoli", "Aeromonas"]
for i=[1,2,3,4,5,6]
    p = scatter(Data[i][1], Data[i][2], yerror=Data[i][3], markersize=2, color=:black, alpha=0.3, title=strain_names[i])
    plot!(sols[i], linewidth=2, color=:red, legend=false, grid=false, xlabel="Time [hr]", ylabel="Height [μm]")
    if i in [1, 2, 5]
        my_box = (1, bbox(0.1, 0.05, 0.4, 0.4))
    else
        my_box = (1, bbox(0.6, 0.55, 0.4, 0.4))
    end
    scatter!(sols[i](Data[i][1]), Data[i][2], label=false, markerstrokecolor=:red, 
            markersize=1, color=:red, grid=false, subplot=2,
            inset = my_box)
    plot!([0, maximum(Data[i][2])], [0, maximum(Data[i][2])], color=:black, alpha=0.5, 
           subplot=2, label=false, linewidth=1.5, linestyle=:dash,
           xticks=[], yticks=[])
    annotate!((0.65, 0.1), text("Measured", 8), subplot=2)
    annotate!((0.1, 0.77), text("Model", 8, rotation=90), subplot=2)

    append!(my_plots, [p])
end
plot(my_plots[1], my_plots[2], my_plots[3], my_plots[4], my_plots[5], my_plots[6], size=(900, 550), left_margin=3mm, dpi=500)
savefig("figs/figs_temp/fig4.png")