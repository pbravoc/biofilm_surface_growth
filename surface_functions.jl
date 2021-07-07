using ModelingToolkit, DifferentialEquations
using Plots, ColorSchemes, Plots.PlotMeasures
using DataFrames, StatsPlots, Arrow
using StatsBase

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

G(x, xstar) = x < xstar ? x : xstar # if below xstar double, otherwise max speed
@register G(x, xstar)

function interface_limited()
    @parameters t α β ϵ Kc
    @variables b(t) c(t)
    D = Differential(t)
    u0 = [b => 0.7, c => 150.0]
    p  = [α => 1.168, β =>0.063, ϵ => 0.012,
          Kc => 140]
    eqs = [D(b) ~ α .*(c/(Kc+c)) *G(b, 20) - β*b,
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

function z_histogram(tf)
      T = length(tf.Profile)
      img = zeros(T, length(tf.Profile[1]))
      hchange = zeros(T-1, length(tf.Profile[1]))
      for i=1:T-1
      hchange[i,:] = (tf.Profile[i+1] - tf.Profile[i]) / (tf.Time[i+1] - tf.Time[i])
      img[i,:] = tf.Profile[i]
      end
      img[T,:] = tf.Profile[T]
      l, r = findall(x->!isnan(x), img[1,:])[1], findall(x->!isnan(x), img[1,:])[end]
      x = vec(img[1:T-1, l:r])
      y = vec(hchange[1:T-1, l:r])
      h = StatsBase.fit(Histogram, (x, y), nbins=50)
      return h
end

function interface_deriv(x, p)
      if_shape(x, p) = p[2]*G.(x, p[1]) - p[3]*x
      return if_shape(x, p)
  end
