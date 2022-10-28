using ModelingToolkit, MethodOfLines, DifferentialEquations, DomainSets
using Distributions 
using DataFrames, CSV
using Interpolations, NumericalIntegration

function interpol(x, y, radius)
    interp = LinearInterpolation(x, y, extrapolation_bc=Line()) 
    x_inter = Array(0:0.01:radius)
    y_inter = [interp(x) for x in x_inter]
    return x_inter, y_inter
end

function diffusion_monod(x_max, t_max, D_biofilm, f, k)
    @parameters x t
    @variables c(..)
    Dt = Differential(t)        # Time derivative
    Dx = Differential(x)        # First spatial der
    Dxx = Differential(x)^2     # Second spatial der

    t_min = 0.0                 # Starting point
    x_min = 0.0                 # μm
    #x_max = R                  # μm
    N = 201                     # discretization 'bins', ODD NUMBER
    dx = (x_max-x_min)/N        # bin width
    c_outside = 1.0
    eqs = [Dt(c(x,t)) ~ D_biofilm*Dxx(c(x, t)) - f*c(x, t)/(k+c(x, t))]

    domains = [x ∈ Interval(x_min, x_max),
               t ∈ Interval(t_min, t_max)]

    bcs = [c(x, 0) ~ c_outside,             # Start at 0
           c(x_min, t) ~ c_outside,                   # Nutrient source     
           Dx(c(x_max, t)) ~ 0.0]               # No flow
    @named pdesys = PDESystem(eqs, bcs, domains,[x,t],[c(x,t)])
    discretization = MOLFiniteDifference([x => dx], t)
    prob = discretize(pdesys, discretization)
    sol = solve(prob, Rodas5(), saveat=1.0)
    x_array = Array(get_discrete(pdesys, discretization)[x])[2:end-1]
    y_array = sol.u[end]
    interp = LinearInterpolation(x_array, y_array, extrapolation_bc=Line()) 
    x_inter = Array(0:0.01:x_max)
    y_inter = interp(x_inter)
    y_monod = y_inter ./ (y_inter .+ k)
    growth_term = integrate(x_inter, y_monod)
    return x_inter, y_inter, y_monod, growth_term
end

function approx_func(x, L)
    if x < L 
        return x 
    else 
        return L 
    end
end

D = 25.0
f = 0.2
k = 0.01
my_sol = diffusion_monod(50.0, 150.0, D, f, k)
L = sqrt(D*k/f)

##
sum_growth = [integrate(my_sol[1][1:i], my_sol[3][1:i]) for i=1:length(my_sol[1])]
my_approx = [approx_func(x, 16.75) for x in my_sol[1]]

df = DataFrame("X"=>my_sol[1], "Concentration"=>my_sol[2]./maximum(my_sol[2]), 
               "Monod"=>my_sol[3] ./ maximum(my_sol[3]), 
               "Sum_growth"=>sum_growth ./ maximum(sum_growth),
               "Approximation"=>my_approx ./ maximum(my_approx))
using Plots, StatsPlots, ColorSchemes
@df df plot(:X, :Concentration, color=:gray,linestyle=:dash, linewidth=2.5,
                 label="Concentration")
@df df plot!(:X, :Monod, color=:gray, linewidth=2.5, label="Monod")
@df df plot!(:X, :Sum_growth, color=ColorSchemes.okabe_ito[1], linewidth=2.5, 
             label="Total Growth")
@df df plot!(:X, :Approximation, color=ColorSchemes.okabe_ito[1], 
             linewidth=2.5, linestyle=:dash, label="Approximation")
plot!(xticks=([0, 10, 20, 30, 40, 50], ["0L", "L", "2L", "3L", "4L", "5L"]), 
      xlim=(0, 40), grid=false, legend=:right, xlabel="Distance from interface",
      ylabel="Normalized value")
#CSV.write("data/sims/monod_diffusion.csv", df)
#@df df plot(:X, :Approximation - :Sum_growth)
