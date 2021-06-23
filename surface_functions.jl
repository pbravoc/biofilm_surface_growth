using ModelingToolkit, DifferentialEquations
using Plots, ColorSchemes, Plots.PlotMeasures
using DataFrames, StatsPlots, Arrow

function nutrient_limited()
    @parameters t α β Kb Kc
    @variables b(t) c(t)
    D = Differential(t)

    u0 = [b => 0.01, c => 1.0]
    p  = [α => 36.0, β => 24/3,
          Kb => 1.0, Kc => 3]
    eqs = [D(b) ~ α*b*(c/(Kc+c))*(1-b/Kb),
          D(c) ~ -β*b*(c/(Kc+c))]
    tspan = (0.0, 2.0)

    prob = ODEProblem(ODESystem(eqs),u0,tspan,p)
    sol = solve(prob,Tsit5())
    return sol 
end

function interface_limited()
    @parameters t α β ϵ Kc
    @variables b(t) c(t)
    D = Differential(t)

    u0 = [b => 0.7, c => 150.0]
    p  = [α => 1.168, β =>0.063, ϵ => 0.012,
          Kc => 140]
    eqs = [D(b) ~ α .*(c/(Kc+c)) *G(b) - β*b,
          D(c) ~ -α*ϵ*b*(c/(Kc+c))]

    prob = ODEProblem(ODESystem(eqs),u0,tspan,p)
    sol = solve(prob,Tsit5(), saveat=0.1)
    return sol 
end

function experimental_heights(strain, replicate)
    df = DataFrame(Arrow.Table("/home/pablo/Biofilms/Data/radialv2.arrow"));
    tf = df[(df.Strain .== strain) .& (df.Replicate .<= replicate) , :]
    center_ids = Int.((df.homelandL+df.homelandR)/2)
    df.height = [df.Profile[i][center_ids[i]] for i=1:size(df)[1]]
    t_experiment, h_experiment = tf.Time, tf.Z
    return t_experiment, h_experiment
end

G(x) = x < 20 ? x : 20 # if below xstar double, otherwise max speed
@register G(x)