##
using Plots, SpecialFunctions

function c(x, t)
    L = sqrt(D*t)
    return erfc.(x / L)
end

function int_c(h)
    L = sqrt(D*dt)
    t1 = L*(-exp.(-h^2 / (L^2)))
    t2 = sqrt(π) * h * erfc(h/L)
    return (t1 + t2 + L)
end

x = Array(0.0:0.05:3.0)
D = 0.3
dt = 1.0
y = exp.(-D*x)
myc = c(x, dt)
a = int_c(0.1)
sol_dis = cumsum(c(x, dt)) .- 1
sol_ana = int_c.(x)
plot(x, myc, xlabel="Distance", grid=false, label="Concentration")
plot!(x, sol_dis ./ maximum(sol_dis), label="Cumulative (sum)")
plot!(x, sol_ana ./ maximum(sol_ana), label="Cumulative (analyitical)")
vline!([sqrt(D*dt)], color=:black, linestyle=:dash, legend=:right, label="Diffusion length")

## Two interfaces
function c(x, D, t)
    L = sqrt(D*t)
    return erfc.(x / L)
end

x = Array(0.0:0.05:7.0)

dt = 1.0
D_oxy = 20.0
D_glu = 6.70
c_oxy = reverse(c(x, D_oxy, dt))
c_glu = c(x, D_glu, dt)
min_c = [minimum([c_oxy[i], c_glu[i]]) for i=1:length(c_oxy)]
mul_c = c_oxy .* c_glu
p1 = plot(x, [c_oxy, c_glu], label=["Oxygen" "Glucose"], legend=:top, grid=false, xlabel="Distance from substrate", ylabel="Concentration")
p2 = plot(x, [mul_c / maximum(mul_c), min_c / maximum(min_c)], label=["Multiplicative" "Minimum"], c=[3 4], grid=false, xlabel="Distance from substrate", ylabel="Joint concentration")
plot(p1, p2, layout=(2,1), dpi=300)
savefig("figs/figs_temp/spatial_concentrations.svg")
##
t_oxy = 1.03*300*300/D_oxy 
t_glu= 1.03*300*300/D_glu
t_oxy