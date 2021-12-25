# This is the code to fit the experimental data 
# to our interface-geometrical model 

using Plots, DataFrames, CSV
using Turing, DifferentialEquations 

# Load up a trajectory of data 
df = DataFrame(CSV.File("biofilm_surface_growth/data/timelapses/database.csv"))

strain, repli = "SN503", "A"
tf = filter(row -> row.Replicate .== "A" && row.Strain .== strain, df)
scatter(tf.Time, abs.(tf.mid_height), color=:black, grid=false, label=string(strain,repli), 
        xlabel="Time (hr)", ylabel="Height (μm)", legend=:topleft)

## Now we define the interface model
G(z, zstar) = z .< zstar ? z : zstar 
@register G(z, zstar)

"""
Interface limited model for biofilm growth,
assumes that nutrients are infinite
"""
function interface_limited(du, u, p, t)
    h = u[1] 
    α, β, hstar = p
    du[1] = α .*G.(h, hstar) .- β .* h 
    return du
end

p = [0.96, 0.057, 12.0]
u0 = 0.12
prob = ODEProblem(interface_limited, [u0], (0.0, 50.0), p) # Set the problem
sol = solve(prob, saveat=0.01)

plot!(sol, color=1, linestyle=:dash, label="Initial guess")
## Now with probabilities

##
Turing.setadbackend(:forwarddiff)

@model function fit_IM(data, prob)
    σ ~ InverseGamma(2, 3) # ~ is the tilde character
    α ~ truncated(Normal(1.5,0.5),0.5,2.5)
    β ~ truncated(Normal(1.2,0.5),0,2)
    γ ~ truncated(Normal(3.0,0.5),1,4)
    δ ~ truncated(Normal(1.0,0.5),0,2)

    p = [α,β,γ,δ]
    prob = remake(prob, p=p)
    predicted = solve(prob,Tsit5(),saveat=0.1)

    for i = 1:length(predicted)
        data[:,i] ~ MvNormal(predicted[i], σ)
    end
end

model = fitlv(odedata, prob1)

##
α = Truncated(Normal(p[1],0.2), p[1]-0.5,p[1]+0.5)
β = Truncated(Normal(p[2],0.005), p[2]-0.01,p[2]+0.01)
h_star = Truncated(Normal(p[3],2), p[3]-10,p[3]+10)
p_prob = [α,β,h_star]

##
print("remaking the problem lol")
prob_prob = remake(prob, p_prob)
print("remade the problem lol")
##
prob = ODEProblem(interface_limited, [u0], (0.0, 50.0), p_prob) # Set the problem
##
predicted = solve(prob, Tsit5(), saveat=0.5)

##

0.1 < h_star 