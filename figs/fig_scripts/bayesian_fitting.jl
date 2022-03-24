using Turing, Distributions, DifferentialEquations
using MCMCChains, Plots, StatsPlots
using Random, SpecialFunctions
using DataFrames, CSV
Random.seed!(14);

function G(h, L)
    t1 = L*(-exp.(-h^2 / (L^2)))
    t2 = sqrt(π) * h * erfc(h/L)
    return (t1 + t2 + L)
end
function interface(du, u, p, t)
    h = u[1] 
    α, β, L = p
    du[1] = α*G(h, L) - β*h 
    return du
end

@model function interface_model(data, prob1, p_guess, save_times)
    α ~ Normal(p_guess[1], p_guess[1]/4) # Growth 
    β ~ Normal(p_guess[2], p_guess[2]/4) # Decay
    #L ~ Normal(15.0, 1.0) # Diffusion Length
    L = p_guess[3]                       # Fixed doesnt break!
    σ ~ InverseGamma(10,10)              # Noise

    p = [α, β, L]
    prob = remake(prob1, p=p)            # Update from sample
    predicted = solve(prob, saveat=save_times, save_idxs=1)
    
    for i = 1:length(predicted)          # Populate array
        data[i] ~ Normal(predicted[i], σ) # predicted[i][2] is the data for y - a scalar, so we use Normal instead of MvNormal
    end
end

## Fake data
t_sample = Array(0.5:1.5:48.0)
p = [0.8, 0.1, 15.0]
u0 = [0.2]
prob1 = ODEProblem(interface,u0,(0.0,50.0),p)
sol1 = solve(prob1,Tsit5(),saveat=t_sample, save_idxs=1)
odedata = abs.(Array(sol1) + 1.0 * randn(size(Array(sol1))))
plot(sol1, alpha = 0.3, linewidth=2, legend = false); 
scatter!(sol1.t, odedata)
##
model = interface_model(sol1.u, prob1, p, t_sample)
chain = sample(model, NUTS(0.65), MCMCThreads(), 1000, 3)
plot(chain)

## Now with real data!
Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.strain .== "bgt127" && 
               x.time .<= 48 && x.replicate in ["A", "B", "C"], Df)
@df df scatter(:time, :avg_height, yerror=:std_height, grouping=:replicate)
p = [0.7590526257703044, 0.05155605894558507, 15.243745168270465]
idxs = sortperm(df.time)
t_data, h_data = df.time[idxs], df.avg_height[idxs]
scatter(t_data, h_data)
u0 = [0.01]
prob1 = ODEProblem(interface,u0,(0.0,50.0),p)
sol = solve(prob1, saveat=0.5)
plot!(sol, color=:black, linewidth=2)
##
model = interface_model(h_data, prob1, p, t_data)
chain = sample(model, NUTS(0.65), MCMCThreads(), 1000, 3)

##
alpha_values = vec(Array(chain[:,:α,:]))
beta_values = vec(Array(chain[:,:β,:]))
α = quantile(alpha_values, [0.025, 0.5, 0.975])
β = quantile(beta_values, [0.025, 0.5, 0.975])

scatter(t_data, h_data)
for i=1:3
    prob = remake(prob1, p=[α[i], β[i], 15.243])
    sol = solve(prob)
    plot!(sol, linewidth=2)
end
print(α ./ β * p[3])
plot!(xlabel="Time")
##
plot(chain, size=(800, 600), left_margin=3mm)
savefig("figs/bayesian_fitting/chain_plot.svg")
##
corner(chain, size=(800, 600), left_margin=3mm)
savefig("figs/bayesian_fitting/chain_corner.svg")
