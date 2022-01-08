using DifferentialEquations, DiffEqFlux, Plots
using DataFrames, CSV, SpecialFunctions

G(z, zstar) = z < zstar ? z : zstar 

function diff_nutrient(h, L)
    t1 = L*(-exp.(-h^2 / (L^2)))
    t2 = sqrt(π) * h * erfc(h/L)
    return (t1 + t2 + L)
end

function interface_limited(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

function diffusion_limited(du, u , p, t)
    h = u[1]
    α, β, L = p 
    du[1] = α*diff_nutrient(h, L) - β*h 
    return du
end

function if_fit(time_data, x_data, p)
    u0 = x_data[1]  
    prob = ODEProblem(interface_limited, [u0], (0.0, 50.0), p) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=time_data) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- x_data)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, p,
                                        maxiters=100)
    return result_ode.u
end

function di_fit(time_data, x_data, p)
    u0 = x_data[1]  
    prob = ODEProblem(diffusion_limited, [u0], (0.0, 50.0), p) # Set the problem
    function loss(p)
        sol = solve(prob, Tsit5(), p=p, saveat=time_data) # Force time savings to match data
        sol_array = reduce(vcat, sol.u)
        loss = sum(abs2, sol_array .- x_data)
        return loss, sol
    end
    result_ode = DiffEqFlux.sciml_train(loss, p,
                                        maxiters=100)
    return result_ode.u
end

df = DataFrame(CSV.File("data/timelapses/database.csv"))
my_strain, my_replicate = "BGT127", "A"
tf =  filter(row -> row.Replicate .== my_replicate && 
             row.Strain .== my_strain, df);
p = [1.0, 0.05, 15.0]

p1 = scatter(tf.Time, tf.loess_height, label="Experimental Data", 
        grid=false, title=string(my_strain, my_replicate))
##
prob = ODEProblem(diffusion_limited, [tf.loess_height[1]], (0.0, 50.0), p) # Set the problem
sol = solve(prob, saveat=tf.Time)

##
iffit = if_fit(tf.Time, tf.loess_height, p)
difit = di_fit(tf.Time, tf.loess_height, p)

##
probdi = ODEProblem(diffusion_limited, [tf.loess_height[1]], (0.0, 50.0), difit) # Set the problem
soldi = solve(probdi, saveat=tf.Time)
probif = ODEProblem(interface_limited, [tf.loess_height[1]], (0.0, 50.0), iffit) # Set the problem
solif = solve(probif, saveat=tf.Time)
err_di = sqrt(mean(abs2, reduce(vcat, soldi.u) .- tf.avg_height))
err_if = sqrt(mean(abs2, reduce(vcat, solif.u) .- tf.avg_height))

#h_pred = round(mlfit[1]*mlfit[3]/mlfit[2], digits=1)
scatter(tf.Time, tf.avg_height, label="Experimental Data", legend=:topleft, 
        grid=false, title=string(my_strain, my_replicate), color=:black)
plot!(soldi, color=1, linewidth=2, xlabel="Time (hr)", label=string("Diffusion, err=", round(err_di, digits=2)))
plot!(solif, color=2, linewidth=2, xlabel="Time (hr)", label=string("Interface, err=", round(err_if, digits=2)))

##
soldi = solve(probdi, saveat=0.1)
hdi = reduce(vcat, soldi.u)
dh_di = (hdi[2:end]-hdi[1:end-1])/0.1
solif = solve(probif, saveat=0.1)
hif = reduce(vcat, solif.u)
dh_if = (hif[2:end]-hif[1:end-1])/0.1

##
p1 = scatter(tf.Time, tf.local_slope, color=:black, alpha=0.7)
p1 = plot!(soldi.t[1:end-1], dh_di, color=1, linewidth=2)
p1 = plot!(solif.t[1:end-1], dh_if, color=2, linewidth=2)

p2 = scatter(tf.loess_height, tf.local_slope, color=:black, alpha=0.7)
p2 = plot!(hdi[1:end-1], dh_di, color=1, linewidth=2)
p2 = plot!(hif[1:end-1], dh_if, color=2, linewidth=2)

plot(p1, p2, size=(700, 300))