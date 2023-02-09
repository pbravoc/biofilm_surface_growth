using ModelingToolkit, MethodOfLines, DifferentialEquations, DomainSets
using Distributions 
using DataFrames, CSV
using Interpolations, NumericalIntegration
using Plots, StatsPlots, ColorSchemes
using Optim 

function interpol(x, y, radius)
    interp = LinearInterpolation(x, y, extrapolation_bc=Line()) 
    x_inter = Array(0:0.01:radius)
    y_inter = [interp(x) for x in x_inter]
    return x_inter, y_inter
end

function diffusion_monod(x_max, t_max, D_biofilm, λ, c_half, c_outside)
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
    eqs = [Dt(c(x,t)) ~ D_biofilm*Dxx(c(x, t)) - λ*c(x, t)/(c_half+c(x, t))]

    domains = [x ∈ Interval(x_min, x_max),
               t ∈ Interval(t_min, t_max)]

    bcs = [c(x, 0) ~ c_outside,             # Start at 0
           c(x_min, t) ~ c_outside,                   # Nutrient source     
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
    y_monod = y_inter ./ (y_inter .+ c_half)
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

# oxygen parameters
D = 25.0                    # μm^2 /s
λ = 1.1e2                 # μM / s
k = 1.0                # μM
c_0= 250.0            # μM
dimensionless = sqrt(D*k / λ)
dimensionless2 = sqrt(c_0*D*k / λ)
##
my_sol = diffusion_monod(30.0, 1e5, D, λ, k, c_0)
sum_growth = [integrate(my_sol[1][1:i], my_sol[3][1:i]) for i=1:length(my_sol[1])]
f(L) = sqrt(mean((sum_growth./maximum(sum_growth) .- approx_func.(my_sol[1], L[1])).^2))
L_best = Optim.minimizer(Optim.optimize(f, [0.0], [100.0], [30.0]))[1]
my_approx = [approx_func(x, L_best) for x in my_sol[1]]
df = DataFrame("X"=>my_sol[1], "Concentration"=>my_sol[2]./maximum(my_sol[2]), 
               "Monod"=>my_sol[3] ./ maximum(my_sol[3]), 
               "Sum_growth"=>sum_growth ./ maximum(sum_growth),
               "Approximation"=>my_approx ./ maximum(my_approx))
@df df plot(:X, :Concentration, color=:gray,linestyle=:dash, linewidth=2.5,
                 label="Concentration")
@df df plot!(:X, :Monod, color=:gray, linewidth=2.5, label="Monod")
@df df plot!(:X, :Sum_growth, color=ColorSchemes.okabe_ito[1], linewidth=2.5, 
             label="Total Growth")
@df df plot!(:X, :Approximation, color=ColorSchemes.okabe_ito[1], 
             linewidth=2.5, linestyle=:dash, label="Approximation")
vline!([dimensionless2], color=:red, linewidth=2)
#plot!(xticks=([0.5*L_best, L_best, 1.5*L_best, 2*L_best], ["0.5L", "L", "1.5L", "2L"]))
plot!( grid=false, legend=:right, xlabel="Distance from interface",
      ylabel="Normalized value")
##
#CSV.write("data/sims/monod_diffusion.csv", df)
#@df df plot(:X, :Approximation - :Sum_growth)
## Behavior over different c halves 
D = 800
λ = 1.3e3
c_outside = 100.0
L = sqrt(D*c_half / (λ*ρ))
c_half_list = Array(5:5:100)
L_list = zeros(length(c_half_list))
error_list = zeros(length(c_half_list))
SG_list = []
MONOD_list = []
for i in 1:length(c_half_list)
    my_sol = diffusion_monod(60.0, 1e5, D, λ, c_half_list[i], c_outside)
    sum_growth = [integrate(my_sol[1][1:i], my_sol[3][1:i]) for i=1:length(my_sol[1])]
    res = Optim.optimize(f, [0.0], [100.0], [30.0])
    L_list[i] = res.minimizer[1]
    error_list[i] = res.minimum
    append!(SG_list, [sum_growth ./ maximum(sum_growth)])
    append!(MONOD_list, [my_sol[3]/maximum(my_sol[3])])
end
##
res = Optim.optimize(f, [0.0], [100.0], [30.0])
##
p1 = plot(c_half_list, error_list, ylim=(0, 0.06), legend=false, c=:black)
plot!(xlabel="Monod Constant [μM]", ylabel="RMSE", right_margin=13mm,
      xticks=[10, 30, 50, 70, 90])
##
p2 = plot(my_sol[1], MONOD_list, line_z = transpose(c_half_list), legend=false, colorbar=true,
          c=:viridis, xlim=(0,40), xlabel="X [μm]", ylabel="Growth", colorbartitle="Monod c.")
p3 = plot(my_sol[1], SG_list, line_z = transpose(c_half_list), legend=false, colorbar=true,
          c=:viridis, xlim=(0,40), xlabel="X [μm]", ylabel="Total Growth", colorbartitle="Monod c.")
##
l = @layout([a{0.4w} b{0.4w} c{0.2w}])
plot(p2, p3, p1, layout=l, size=(700, 150), grid=false, dpi=500, bottom_margin=5mm, left_margin=3mm)
savefig("figs/monod_constant_approximation.svg")
