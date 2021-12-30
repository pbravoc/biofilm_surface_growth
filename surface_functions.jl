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

function interface_deriv(x, p)
      if_shape(x, p) = p[2]*G.(x, p[1]) - p[3]*x
      return if_shape(x, p)
  end

  " Substract a polynomial of degree *deg* from the array y"
  function subpoly(y, deg)
      x = Array(range(1, length(y), step=1))  # temp x array for fitting
      flat = fit(x, y, deg)                   # fitting 
      myfit = flat.(x)                        # lower dimension
      sub = y-myfit                           # substract 
      return sub
  end
  
  "Returns frequency and FFT of a pair of arrays *x* and *y*"
  function fourer_decomp(x, y)
      xs = x[2]-x[1]                          # get distance
      F = fft(y) |> fftshift                  # get amplitudes 
      freqs = fftfreq(length(x), 1/xs) |> fftshift    # get frequencies
      return freqs, F
  end
  
  "Returns a mirror of *v* displaced by *s*, NaNs on out of range"
  function mirror(v, s)
      mirr_vector = circshift(v, s)           # PBC shift the vector
      if s>=0                                 # if change is positive
          mirr_vector[1:s] .= NaN             # NaN beginning
      else                                    # if its negative
          mirr_vector[end+s+1:end] .= NaN     # NaN end
      end
      return mirr_vector
  end
  
  """
  Updates first *M1* and second moments *M2* of the 
  heights array *v* across a range of distances *r*
  """
  function update_moments!(v, r, C, M1, M2)
      for r_current in (r, -r)
          vect = mirror(v, r_current)
          idx = isnan.(vect)                  # values that are NaNs
          C .+= .! idx                        # Only count values with data
          vect[idx] .= 0                      # NaNs->0 to not remove data
          M1 .+= vect                         # Update first moment
          M2 .+= vect .^ 2                    # Update second moment
      end
      return C, M1, M2
  end
  
  """
  Calculates the local standard deviation of the
  array *v* on a set of distances *R*
  """
  function lstd(v, R)
      C = zeros(length(v))                    # empty counts 
      M1 = zeros(length(v))                   # empty moment 1
      M2 = zeros(length(v))                   # empty moment 2
      LSD = zeros(length(R)-1)                # local standard deviation
      for i in 1:(length(R)-1)                # loop for all the distances
          C, M1, M2 = update_moments!(v, R[i+1], C, M1, M2)   # update
          ex1 = M1 ./ C                       # Real moment 1 by normalizing
          ex2 = M2 ./ C                       # Real moment 2 by normalizing
          inside = ex2 .- (ex1 .^ 2)          # to make other operation shorter
          LSD[i] = NaNMath.mean(sqrt.(abs.(inside))) # standard deviation
      end                                     # abs solves negative zeroes from comp
      return LSD
  end
  
  """
  Box-counting algorithm for the array *y*, with a spatial
  resolution *R* and at sizes *S*
  """
  function get_counts(y,rx, S)
      x = Array(1:length(y))
      R = [1, rx]   # Resolution
      M = [length(y), NaNMath.maximum(y)]  # Maximums
      m = [0, NaNMath.minimum(y)]          # Minimums
      rs = [R*s for s in S]
      bn = Array(reduce(hcat, [Int.((M-m) .÷ r) for r in rs])')
      bn[bn .== 0] .= 1
      bs = Array(reduce(hcat, [(M-m)./bn[i,:] for i in 1:size(bn)[1]])')
      counts = zeros(length(S))
      for i in 1:length(S)
          c = bn[i,2]
          C = Int.(c*(x .÷ bs[i,1]) + (y .÷ bs[i,2]))
          idx = .! isnan.(C)
          counts[i] = length(unique(C[idx]))
      end
  return counts
  end
  
  @. pline(x,p) = p[1] + p[2] * x
  function plinefit(x,y)
      myfit = curve_fit(pline, x, y, [1.0, 1.0])
      return myfit.param
  end
  
  function log_weights(x)
      xarray = log10.(x)
      xstart, xend = xarray[1:end-1], xarray[2:end]
      xdif = xend-xstart
      return push!(xdif, xdif[end])
  end
  
  function fit_rough(x, y)
      x, y = log10.(x)[8:400], log10.(y)[8:400]
      we = log_weights(x)
      we = we./sum(we)
      return curve_fit(pline, x, y, we, [-3.0, 1.0]).param[2]
  end

  function d_height(df)
    h_change = zeros(size(df)[1]) 
    h_change .= NaN
    if size(df)[1] > 2
        h_change[1:end-1] = (df.avg_height[2:end]-df.avg_height[1:end-1]) ./ 
                            (df.Time[2:end]-df.Time[1:end-1])
    end
    return h_change
end

function smooth_heights(df, time_window)
    loess_height, slope = zeros(size(df)[1]), zeros(size(df)[1])
    loess_height .= NaN
    slope .= NaN
    if size(df)[1] > 2
        per_window = time_window / (df.Time[end]-df.Time[1])
        model = loess(df.Time, df.avg_height, span=per_window, degree=1)
        loess_height = [predict(model, Float64(df.Time[i])) for i=1:length(df.Time)]
        slope = zeros(length(df.Time))
        dt = 0.1
        t = df.Time[1]
        slope[1] = (predict(model, Float64(t+dt))-predict(model,Float64(t))) / dt
        t = df.Time[end]
        slope[end] = (predict(model, Float64(t))-predict(model,Float64(t-dt))) / dt
        for i=2:length(df.Time)-1
            t = df.Time[i]
            slope[i] = (predict(model, Float64(t+dt))-predict(model,Float64(t-dt))) / (2*dt)
        end
    end
    return loess_height, slope
end