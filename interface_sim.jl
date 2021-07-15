using DifferentialEquations, DiffEqOperators
using CSV, DataFrames
using Plots, Plots.Measures

G(z, zstar) = z < zstar ? z : zstar 
@register G(z, zstar)

"""
Spatial interface colony simulation
"""
function biofilm(u0, x, pars, tspan)
    α, β, ν, ucrit= pars
    Δ = CenteredDifference(2, 2, x[2], length(x))
    bc = PeriodicBC(Float64)
    f(u,p,t) =  ν*Δ*bc*u + α*G.(u, ucrit) - β*u
    return ODEProblem(f, u0, (tspan[1], tspan[end]))
end


f = DataFrame(CSV.File("data/sims/merger_initial.csv", header=false) )
x = f[:,1]
z = abs.(f[:,2])
z[y .< 0.065] .= 0
pars = [1.0, 0.069, 60.0, 12.0]            
t = Array(0:1.0:20)                                    # Time range
prob = biofilm(z, x, pars, t)
sol = solve(prob, saveat=0.1);

anim = @animate for i=1:100
    plot(x, sol.u[i], ylim=(0, 80), color=:black, label=false, xlabel="X (μm)", ylabel="Height (μm)", dpi=200)
    scatter!([x[750]], [sol.u[i][750]], color=1, label=false)
    scatter!([x[1710]], [sol.u[i][1710]], color=2, label=false)
    scatter!([x[1210]], [sol.u[i][1210]], color=3, label=false, title=string.("t = ", sol.t[i], "hr"))
end

gif(anim, "figs/sims/merger.gif", fps=15)

l_peak = [sol.u[i][750] for i =1:100]
r_peak = [sol.u[i][1710] for i =1:100]
m_peak = [sol.u[i][1210] for i =1:100]

get_dif(y,t) = (y[2:end] - y[1:end-1]) ./ (t[2:end] - t[1:end-1])
p1 = plot(sol.t[1:100], [l_peak, r_peak, m_peak], labels=["Left" "Right" "Merger"], marker=:circle, legend=:topleft,
     xlabel="Time (hr)", ylabel="Height (μm)")
p2 = plot(l_peak[1:99], get_dif(l_peak, sol.t[1:100]), marker=:circle, label="Left", xlabel="Height (μm)", ylabel="Height change (μm/hr)")
plot!(r_peak[1:99], get_dif(r_peak, sol.t[1:100]), marker=:circle, label="Right")
plot!(m_peak[1:99], get_dif(m_peak, sol.t[1:100]), marker=:circle, label="Merger")

plot(p1, p2, size=(600 ,300), dpi=200)
savefig("figs/sims/merger_tracked_hegihts.svg")