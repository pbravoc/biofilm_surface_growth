using ModelingToolkit, MethodOfLines, DifferentialEquations, DomainSets
using Distributions 
using DataFrames, CSV
using Interpolations, NumericalIntegration
using Plots, StatsPlots, ColorSchemes
using Optim 

function diffusion_monod(x_max, t_max, pars)
    @parameters x t
    @variables c(..)
    Dt = Differential(t)        # Time derivative
    Dx = Differential(x)        # First spatial der
    Dxx = Differential(x)^2     # Second spatial der

    t_min = 0.0                 # Starting point
    x_min = 0.0                 # μm
    N = 201                     # discretization 'bins', ODD NUMBER
    dx = (x_max-x_min)/N        # bin width
    eqs = [Dt(c(x,t)) ~ pars["D"]*Dxx(c(x, t)) - pars["λ"]*c(x, t)/(pars["k"]+c(x, t))]

    domains = [x ∈ Interval(x_min, x_max),
               t ∈ Interval(t_min, t_max)]

    bcs = [c(x, 0) ~ pars["c₀"],             # Start at 0
           c(x_min, t) ~ pars["c₀"],                   # Nutrient source     
           Dx(c(x_max, t)) ~ 0.0]               # No flow
    @named pdesys = PDESystem(eqs, bcs, domains,[x,t],[c(x,t)])
    discretization = MOLFiniteDifference([x => dx], t)
    prob = discretize(pdesys, discretization)
    sol = solve(prob, Rodas5(), save_on=false, save_start=false, save_end=true)
    x_array = vec(sol[x])
    y_array = vec(sol[c(x,t)])
    interp = LinearInterpolation(x_array, y_array) 
    x_inter = Array(0:0.01:x_max)
    y_inter = interp(x_inter)
    y_monod = y_inter ./ (y_inter .+ pars["k"])
    growth_term = integrate(x_inter, y_monod)
    return x_inter, y_inter, y_monod, growth_term
end

function approx_func(x, L)
    if x < L 
        return x / L
    else 
        return L / L 
    end
end

G(h, L) = h < L ? h : L 
function interface(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

function interface2(du, u, p, t)
    h = u[1] 
    α, β = p
    du[1] = α*G.(h, L_best) - β*h 
    return du
end

function interface_monod(du, u, p, t)
    h = u[1]
    α, β = p 
    du[1] = α*M(h) - β*h 
    return du 
end

function loss(p)
    sol = solve(prob, Tsit5(), p=p, saveat=t_exp, save_idxs=1) # Only save height
    loss = sum(abs2, h_exp - sol.u)
    return loss
end 

## Simulate real growth curve from Monod 
parameters = Dict("ρ"=> 1.0, "D"=> 800.0, "λ"=> 1.3e3, "k"=>38, "c₀"=>100.0)
my_sol = diffusion_monod(40.0, 1e5, parameters)
M = LinearInterpolation(my_sol[1], cumsum(my_sol[3]) ./ maximum(cumsum(my_sol[3])))
f(L) = sqrt(mean((cumsum(my_sol[3])./maximum(cumsum(my_sol[3])) .- approx_func.(my_sol[1], L[1])).^2))
L_best = Optim.minimizer(Optim.optimize(f, [0.0], [100.0], [30.0]))[1]
sum_growth = [integrate(my_sol[1][1:i], my_sol[3][1:i]) for i=1:length(my_sol[1])]
my_approx = [approx_func(x, L_best) for x in my_sol[1]]
##
df = DataFrame("X"=>my_sol[1], "Concentration"=>my_sol[2]./maximum(my_sol[2]), 
               "Monod"=>my_sol[3] ./ maximum(my_sol[3]), 
               "Sum_growth"=>sum_growth ./ maximum(sum_growth),
               "Approximation"=>my_approx ./ maximum(my_approx))
p1 = @df df plot(:X, :Concentration, color=:gray,linestyle=:dash, linewidth=2.5,
                 label="Concentration")
@df df plot!(:X, :Monod, color=:gray, linewidth=2.5, label="Monod")
@df df plot!(:X, :Sum_growth, color=ColorSchemes.okabe_ito[2], linewidth=2.5, 
             label="Total Growth")
@df df plot!(:X, :Approximation, color=ColorSchemes.okabe_ito[1], 
             linewidth=2.5, label="Approximation")
plot!(xlim=(0, 40), size=(300, 250))
plot!( grid=false, legend=:right, xlabel="Distance from interface",
      ylabel="Normalized value", dpi=300)
#savefig("ecoli_lserine.png")

## Load experimental data

Df =  DataFrame(CSV.File("data/timelapses/model_predictions_all.csv"))
Df = filter(x-> x.avg_height .> 0 && 
                x.strain .== "jt305", Df)                   
t_exp, h_exp = Df.time, Df.avg_height

my_sol = diffusion_monod(600.0, 1e5, parameters)
M = LinearInterpolation(my_sol[1], cumsum(my_sol[3]) ./ maximum(cumsum(my_sol[3])))
# Fit model with real Monod to experimental data 
u0 = [0.1]
prob = ODEProblem(interface_monod, u0, (0.0, t_exp[end]), [3.0, 0.01])
params_monod = Optim.optimize(loss, [3.0, 0.01]).minimizer
prob = ODEProblem(interface_monod, u0, (0.0, t_exp[end]), params_monod)
sol_monod = solve(prob, saveat=t_exp, save_idxs=1)
prob = ODEProblem(interface, u0, (0.0, t_exp[end]), [0.8, 0.1, 15.0])
params_approx = Optim.optimize(loss, [0.8, 0.1, 15.0]).minimizer
prob = ODEProblem(interface, u0, (0.0, t_exp[end]), params_approx)
sol_approx = solve(prob, saveat=t_exp, save_idxs=1)
prob = ODEProblem(interface2, u0, (0.0, t_exp[end]), [0.8, 0.1])
params_approx2 = Optim.optimize(loss, [0.8, 0.1]).minimizer
prob = ODEProblem(interface2, u0, (0.0, t_exp[end]), params_approx2)
sol_approx2 = solve(prob, saveat=t_exp, save_idxs=1)
rmse(x, y) = sqrt(mean(sum(abs2, x-y)))
rmse_monod = rmse(h_exp, sol_monod)
rmse_approx = rmse(h_exp, sol_approx)
rmse_approx2 = rmse(h_exp, sol_approx2)
##
myc = [ColorSchemes.okabe_ito[2], ColorSchemes.okabe_ito[1], ColorSchemes.okabe_ito[6]]
p2 = @df Df plot(:time, :avg_height, ribbon=:std_height, color=:gray, label="Data",
                 fillalpha=0.4)
plot!(t_exp, [sol_monod.u, sol_approx.u, sol_approx2.u], linewidth=[3 3 1], 
      color = myc', label=["Monod" "min (3fp)" "min (2fp)"], 
      linestyle=[:solid :solid :dash])
plot!(legend=:bottomright, xlabel="Time [hr]", ylabel="Height [μm]")
p3 = @df Df plot(:time, :avg_height - :avg_height, ribbon=:std_height, 
                 fillalpha=0.2, color=:gray, label=false, linestyle=:dash, alpha=0.8)
plot!(t_exp, [sol_monod.u - h_exp, sol_approx.u- h_exp, sol_approx2.u- h_exp], 
      color = myc', linewidth=[3 3 1], label=["Monod" "min (3fp)" "min (2fp)"], 
      linestyle=[:solid :solid :dash])
plot!(legend=false, xlabel="Time [hr]", ylabel="Residual [μm]")
plot(p1, p2, p3, layout=(1,3), grid=false, size=(800, 220), bottom_margin=4mm, left_margin=4mm, dpi=300)
#savefig("monod_approx.svg")
##

