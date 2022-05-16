#= This code loops over the database, and returns the best fits
for the 48 hour data. One for each timelapse, and one for the aggregated
data.
It also calculated the best fit using 48h data + measurements from 
longtime_data.csv
=#
using DataFrames, CSV, CircularArrays
using Statistics, NaNMath, StatsBase
using DifferentialEquations, DiffEqFlux
using Glob

function block_bootstrap(df, n_blocks, block_size)
    indices = CircularArray(1:size(df)[1]-1)
    sample_start = sample(indices, n_blocks)
    sample_temp = reduce(vcat, 
                    [Array(x:x+block_size) for x in sample_start])
    si = [(x % size(df)[1])+1 for x in sample_temp]
    return df[si,:]  
end
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
function fit_data(t_data, h_data, model, u0=[0.1], pguess=[0.8, 0.05, 15.0])# pguess=[0.8, 0.05, 15.0]) # 
    idxs = sortperm(t_data)  # Sort time indexes
    t_data, h_data = t_data[idxs], h_data[idxs]
    prob = ODEProblem(model, u0, (0.0, t_data[end]), pguess)
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=t_data, save_idxs=1) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- h_data)
        return loss, sol
    end 
    result_ode = DiffEqFlux.sciml_train(loss, pguess)
    return result_ode
end 

##
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0, Df)                    # Smaller than 0 values don't make physical sense
df2 = DataFrame(CSV.File("data/timelapses/longtime_data.csv"))
strain_list = unique(Df.strain)
#model_choice, n_parameters = logistic, 2
model_choice, n_parameters = interface, 3
u_list = Array(0.01:0.01:0.5)
P = []
Strain = []
Fit = []
U = []
## Vary the starting condition
for strain in strain_list 
    println(strain)
    df = filter(x-> x.strain .== strain && x.time .<48 &&
                    x.replicate in ["A", "B", "C"], Df)  
    for u in u_list  
        print(u)          
        fit_params = fit_data(df.time, df.avg_height, model_choice, [u])
        append!(P, [fit_params.u])
        append!(Strain, [strain])
        append!(Fit, ["48h"])
        append!(U, [u])
    end
end
##
# Get the best fits for less than 48h 
for strain in strain_list 
    println(strain)
    df = filter(x-> x.strain .== strain && x.time .<48, Df)              
    fit_params = fit_data(df.time, df.avg_height, model_choice)
    append!(P, [fit_params.u])
    append!(Strain, [strain])
    append!(Fit, ["48h"])
end
##
# Get the best fits for each timelapse
for strain in strain_list 
    println(strain)
    for replicate in ["A", "B", "C"]     
        df = filter(x-> x.strain .== strain && x.time .<48 &&
                        x.replicate .==replicate, Df)
        fit_params = fit_data(df.time, df.avg_height, model_choice)
        append!(P, [fit_params.u])
        append!(Strain, [strain])
        append!(Fit, [replicate])
    end
end
# This is to get the 'long time' best fit.
strain_list = ["bgt127", "gob33", "jt305"]
df = filter(x-> x.strain in strain_list && x.time .< 48 &&
                x.replicate in ["A", "B", "C"], Df)
for strain in strain_list 
    tf = filter(x-> x.strain .== strain, df)
    tf2 = filter(x-> x.strain .== strain, df2)
    t = append!(tf.time, tf2.time)
    h = append!(tf.avg_height, tf2.avg_height)
    p = fit_data(t, h, model_choice)
    append!(P, [p.u])
    append!(Strain, [strain])
    append!(Fit, ["long"])
end
##
pf = hcat(DataFrame("strain"=>Strain, "fit"=>Fit, "u0"=>U),
          DataFrame(Matrix(reduce(hcat, P)'), :auto))
## Save to file
CSV.write("data/timelapses/fit_params_u0.csv", pf)

 
##
using Plots, StatsPlots, ColorSchemes, Colors

my_colors = [ColorSchemes.okabe_ito[8], ColorSchemes.okabe_ito[5],
             ColorSchemes.okabe_ito[4], ColorSchemes.okabe_ito[6]]
pf.h = pf.x1 .* pf.x3 ./ pf.x2
@df pf scatter(:u0, :h, group=:strain, ylabel="Predicted Height [μm]", xlabel="Starting Height [μm]",
               size=(350, 300), ylim=(0, 610), color=[my_colors[2], my_colors[4], my_colors[3]]',
               marker=[:diamond :square :circle],
               label = ["Aeromonas" "Yeast (aa)" "E coli"],
               legend=:bottomright, markersize=3, dpi=500)
savefig("figs/fig3/starting_conditions.png")
