using Turing, Distributions, DifferentialEquations
using MCMCChains, Plots, StatsPlots
using Random
Random.seed!(14);
f(u,p,t) = p[1]*u*(1- u/(p[2]))
p = [1.5, 150.0]
u0 = 1.0
prob1 = ODEProblem(f,u0,(0.0,10.0),p)
sol1 = solve(prob1,Tsit5(),saveat=0.1)
odedata = Array(sol1) + 1.0 * randn(size(Array(sol1)))
plot(sol1, alpha = 0.3, legend = false); 
scatter!(sol1.t, odedata)
##
@model function fitlog(data, prob1)
    σ ~ InverseGamma(10,10) # ~ is the tilde character
    α ~ Normal(1.5, 0.2)
    K_h ~ Normal(150.0, 10.0)

    p = [α,K_h]
    prob = remake(prob1, p=p)
    predicted = solve(prob, saveat=0.1)

    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ) # predicted[i][2] is the data for y - a scalar, so we use Normal instead of MvNormal
    end
end

model = fitlog(sol1.u, prob1)
chain = sample(model, NUTS(0.65), 3_000);
#chain = mapreduce(c -> sample(model, NUTS(.65), 3_000), chainscat, 1:3)

##
plot(chain)

##
plot()
chain_array = Array(chain)
for k in 1:300
    resol = solve(remake(prob1, p=chain_array[rand(1:1500), 2:3]),Tsit5(),saveat=0.1)
    plot!(resol, alpha=0.1, color = :red, legend = false)
end
plot!(sol1)

##
plot(Normal(150.0, 10.0))