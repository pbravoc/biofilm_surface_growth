using ModelingToolkit, DifferentialEquations
using Plots, ColorSchemes, Plots.PlotMeasures
using DataFrames, StatsPlots, Arrow

function run_nutrients()
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

function printme(x)
    print(x)
end
