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

Df = jldopen("data/timelapses/profile_database.jld2", "r")["df"];
Pf = DataFrame(Arrow.Table("/home/pablo/Biofilms/scripts/biofilm_surface_growth/data/sims/fit_params.arrow"));

