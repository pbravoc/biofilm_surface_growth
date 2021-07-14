using DifferentialEquations, DiffEqParamEstim
using Optim
using Plots
using DataFrames, Arrow

G(z, zstar) = z < zstar ? z : zstar 
@register G(z, zstar)

"""
Interface limited model for biofilm growth,
assumes that nutrients are infinite
"""
function interface_limited(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α*G.(h, hstar) - β*h 
    return du
end

"""
Fits the experimental data to the interface 
limited interface model. Returns the best parameters
growth, decay and critical height
"""
function fit_model(prob, tdata, zdata)
    cost_function = build_loss_objective(prob, Tsit5(), L2Loss(tdata, zdata),
                    maxiters=100000, verbose=false)
    result_bfgs = Optim.optimize(cost_function, [0.8, 0.04, 14], Optim.BFGS())
    min = result_bfgs.minimizer
    return min 
end

function get_data(Df, strain, repl)
    df =  filter(row -> row.Replicate .== repl && row.Strain .== strain, Df);
    return df.Time, df.mid_height
end

# Set problem + dummy parameters
u0, p = [0.1], [0.9, 0.1, 20] # Dummy starting conditions
prob = ODEProblem(interface_limited, u0, (0.0, 50.0), p) # Set the problem


# Load data
Df = DataFrame(Arrow.Table("data/timelapses/profile_database.arrow"));
replicates = ["A", "B", "C"]

# Fit data
F = zeros(3,3)
T, Z = [], []
for i = 1:3
    t, z = get_data(Df, "BGT127", replicates[i])
    f = fit_model(prob, t, z)
    F[i,:] = f 
    push!(T, t)
    push!(Z, z)
end

# Simulate with fitting parameters
S = []
for i =1:3
    prob = ODEProblem(interface_limited, u0, (0.0, 400.0), F[i,:]) # Set the problem
    sol = solve(prob)
    push!(S, sol)
end

# Plot prediction
plot()
for i=1:3
    plot!(S[i], color=i, linewidth=2, label=replicates[i], xlabel="Time (hr)", ylabel="Height (μm)", xlim=(0, 350), size=(500, 400), legend=:bottomright)
end
@df filter(row -> any(row.Replicate .== ["G", "H", "I"]), Df) scatter!(:Time, :mid_height, color=:red, alpha=0.7, label="Control")
savefig("figs/timelapses/bgt127/fit_output.svg")

