using DifferentialEquations, DiffEqParamEstim, Optim
using DataFrames, JLD2
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Plots.Measures

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
function fit_model(prob, tdata, zdata, par_guess)
    cost_function = build_loss_objective(prob, Tsit5(), L2Loss(tdata, zdata),
                    maxiters=100000, verbose=false)
    result_bfgs = Optim.optimize(cost_function, par_guess, Optim.BFGS())
    min = result_bfgs.minimizer
    return min 
end

function interface_limited_multiplepars(du, u, p, t)
    α, β, hstar = p[1,:], p[2,:], p[3,:]
    my_values = α.*G.(u, hstar) .- β.*u
    du[1] = my_values[1]
    du[2] = my_values[2]
    #print(du)
    return du
end

Df = jldopen("data/timelapses/profile_database.jld2", "r")["df"];

# Set problem + dummy parameters
u0, p = [0.1], [0.8, 0.1, 0.2] # Dummy starting conditions
prob = ODEProblem(interface_limited, u0, (0.0, 50.0), p) # Set the problem

dist = Array(5000:100:30000)
adist = abs.(dist .- 17500)

df = Df[(Df.Strain .== "SN503") .& (Df.Replicate .<= "A") , :];
x = (Array(1:length(df.Profile[1])) .- length(df.Profile[1])/2) * 0.17362 * 1e-3
p1 = plot(x, df.Profile, line_z = df.Time', color=:turbo, label=false, ylim=(0, 220), xlabel="X (mm)", ylabel="Height (μm)", colorbar_title="Time (hr)");
vline!([x[9000], x[12000], x[23000], x[26000]], color=:black, style=:dash, label=false, xlim=(-3, 3), grid=false,  size=(1000, 400), bottom_margin=3mm, left_margin=4mm, dpi=300)
annotate!(x[6000], 200, "Out")
annotate!(x[10500], 200, "CR")
annotate!(x[17500], 200, "Homeland")
annotate!(x[24500], 200, "CR")
annotate!(x[28000], 200, "Out");
plot!(x[dist], params[1,:] .* params[3,:] ./ params[2,:], color=:black, linewidth=2, label=false)
savefig("/home/pablo/Biofilms/scripts/biofilm_surface_growth/figs/sims/spatial/timelapse_prediction.png")

plot(df.Time, z, line_z = abs.(xd)', c=:viridis, label=false, colorbar_title="Distance to center (mm)", grid=false, size=(1000, 400), bottom_margin=3mm, left_margin=4mm, dpi=300, alpha=0.8, xlabel="Time (hr)", ylabel="Height (μm)", ylim=(0, 210))
savefig("/home/pablo/Biofilms/scripts/biofilm_surface_growth/figs/sims/spatial/timelapse_trajectories.png")

scatter(x[dist], params', label=["α" "β" "h*"])
vline!([x[9000], x[12000], x[23000], x[26000]], color=:black, style=:dash, label=false, xlim=(-3, 3), grid=false,  size=(1000, 400), bottom_margin=3mm, left_margin=4mm, dpi=300)
annotate!(x[6000], 18, "Out",  xlabel="X (mm)", ylabel="Value")
annotate!(x[10500], 18, "CR")
annotate!(x[17500], 18, "Homeland", ylim=(-0.5, 20))
annotate!(x[24500], 18, "CR")
annotate!(x[28000], 18, "Out");
savefig("/home/pablo/Biofilms/scripts/biofilm_surface_growth/figs/sims/spatial/fit_absolute.png")

scatter(x[dist], (params ./ params[:,125])', label=["α" "β" "h*"])
vline!([x[9000], x[12000], x[23000], x[26000]], color=:black, style=:dash, label=false, xlim=(-3, 3), grid=false,  size=(1000, 400), bottom_margin=3mm, left_margin=4mm, dpi=300)
annotate!(x[6000], 1.4, "Out",  xlabel="X (mm)", ylabel="Relative Value")
annotate!(x[10500], 1.4, "CR")
annotate!(x[17500], 1.4, "Homeland", ylim=(0.5, 1.5))
annotate!(x[24500], 1.4, "CR")
annotate!(x[28000], 1.4, "Out");
savefig("/home/pablo/Biofilms/scripts/biofilm_surface_growth/figs/sims/spatial/fit_relative.png")