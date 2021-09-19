using DataFrames, Arrow, JLD2
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Plots.Measures
using DifferentialEquations, DiffEqParamEstim
using Optim

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
    result_bfgs = Optim.optimize(cost_function, [0.3, 0.015, 15.5], Optim.BFGS())
    min = result_bfgs.minimizer
    return min 
end

Df = jldopen("data/timelapses/profile_database.jld2", "r")["df"];
param_frame = jldopen("/home/pablo/Biofilms/scripts/biofilm_surface_growth/figs/fig_scripts/fit_params.jld2", "r")["param_frame"];

EC_A = abs.(filter(row -> row.Replicate .== "A" && row.Strain .== "JT305L", Df).avg_height)
EC_B = abs.(filter(row -> row.Replicate .== "B" && row.Strain .== "JT305L", Df).avg_height)
EC_C = abs.(filter(row -> row.Replicate .== "C" && row.Strain .== "JT305L", Df).avg_height)
EC_C[52] += 30
EC_C[53] -= 12
EC_C[78] += 40
tEC_A = filter(row -> row.Replicate .== "A" && row.Strain .== "JT305L", Df).Time[1:end-1]
tEC_B = filter(row -> row.Replicate .== "B" && row.Strain .== "JT305L", Df).Time[1:end-1]
tEC_C = filter(row -> row.Replicate .== "C" && row.Strain .== "JT305L", Df).Time[1:end-1]

Ae_A = abs.(filter(row -> row.Replicate .== "A" && row.Strain .== "BGT127", Df).avg_height)
Ae_B = abs.(filter(row -> row.Replicate .== "B" && row.Strain .== "BGT127", Df).avg_height)
Ae_C = abs.(filter(row -> row.Replicate .== "C" && row.Strain .== "BGT127", Df).avg_height)
tAe_A = filter(row -> row.Replicate .== "A" && row.Strain .== "BGT127", Df).Time
tAe_B = filter(row -> row.Replicate .== "B" && row.Strain .== "BGT127", Df).Time
tAe_C = filter(row -> row.Replicate .== "C" && row.Strain .== "BGT127", Df).Time

tEC, EC = [tEC_A, tEC_B, tEC_C], [EC_A, EC_B, EC_C]
tAe, Ae = [tAe_A, tAe_B, tAe_C], [Ae_A, Ae_B, Ae_C]

pbase = plot()
for i=1:3
    scatter!(tEC[i], EC[i][1:end-1], c=ColorSchemes.Blues_3[i], alpha=0.7, label=false)
    scatter!(tAe[i], Ae[i], c=ColorSchemes.Oranges_3[i], alpha=0.7, label=false)
end
for i=1:2
    if i==1
        fit_stuff = Array(param_frame[1:3,2:4])
        u0 = 0.1
    else
        fit_stuff = Array(param_frame[4:6,2:4])
        u0 = 0.4
    end
    m, s = mean(fit_stuff, dims=1), std(fit_stuff, dims=1)
    avg_vals = m
    low_vals = abs.([m[1]-s[1], m[2]+ s[2], m[3]-s[3]])
    hig_vals = abs.([m[1]+s[1], m[2]- s[2], m[3]+s[3]])
    prob = ODEProblem(interface_limited, [u0], (0.0, 348.0), avg_vals) # Set the problem
    sol_avg = solve(prob, saveat=0.1)
    prob = ODEProblem(interface_limited, [u0], (0.0, 348.0), low_vals) # Set the problem
    sol_low = solve(prob, saveat=0.1)
    prob = ODEProblem(interface_limited, [u0], (0.0, 348.0), hig_vals) # Set the problem
    sol_hig = solve(prob, saveat=0.1)
    ribbon_low = abs.(reduce(vcat, sol_avg.u)- reduce(vcat, sol_low.u))/2
    ribbon_hig = abs.(reduce(vcat, sol_avg.u)- reduce(vcat, sol_hig.u))/2
    if i==1
        pbase = plot!(sol_avg.t, reduce(vcat, sol_avg.u), ribbon = (ribbon_low, ribbon_hig), color=ColorSchemes.Oranges_3[3], fillalpha=0.2, linewidth=2,
        label="Aeromonas")
    else 
        pbase = plot!(sol_avg.t, reduce(vcat, sol_avg.u), ribbon = (ribbon_low, ribbon_hig), color=ColorSchemes.Blues_3[3], fillalpha=0.2, linewidth=2,
        label="E. coli")
    end
end

pbase = @df filter(row -> row.mid_height .>= 230, Df[Df.Strain .== "JT305L", :]) scatter!(:Time, :mid_height, marker=:circle, alpha=0.7, 
            xlabel="Time (hr)", ylabel="Height (μm)", color=ColorSchemes.Blues_3[3], label=false)

pbase = @df filter(row -> any(row.Replicate .== ["G", "H", "I"]), Df[Df.Strain .== "JT305", :]) scatter!(:Time, :avg_height, marker=:circle, alpha=0.7, 
            xlabel="Time (hr)", ylabel="Height (μm)", color=ColorSchemes.Blues_3[3],label=false)

            p1 = @df filter(row -> any(row.Replicate .== ["G", "H", "I"]), Df[Df.Strain .== "BGT127", :]) scatter!(:Time, :avg_height, marker=:circle, alpha=0.7, 
            xlabel="Time (hr)", ylabel="Height (μm)", color=ColorSchemes.Oranges_3[3],label=false)
#plot!(logistic_ecoli, color=ColorSchemes.Blues_3[3], linestyle=:dash, label=false)
#plot!(logistic_aero, color=ColorSchemes.Oranges_3[3], linestyle=:dash, label=false)
plot!([480, 488], [100, 100], color=:black, label="Interface model")
plot!([480, 488], [100, 100], color=:black, linestyle=:dash, label="Logistic model", xlabel="Time (hours)")

p1 = plot(pbase, xlim=(0, 48.65), ylim=(-1, 210), grid=false, legend=:topleft)
p2 = plot(pbase, xlim=(0, 348.65), ylim=(-1, 360), grid=false, legend=:topleft)
plot(p1, p2, size=(1000, 400), left_margin=3mm, bottom_margin=3mm)


