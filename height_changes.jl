using Plots, StatsPlots, Colors
using DataFrames, CSV
using Loess 

function d_height(df)
    h_change = zeros(size(df)[1]) 
    h_change[end] = NaN
    if size(df)[1] > 1
        h_change[1:end-1] = (df.avg_height[2:end]-df.avg_height[1:end-1]) ./ 
                            (df.Time[2:end]-df.Time[1:end-1])
    end
    return h_change
end

function smooth_heights(df, time_window)
    per_window = time_window / (df.Time[end]-df.Time[1])
    model = loess(df.Time, df.avg_height, span=per_window, degree=1)
    loess_height = [predict(model, t) for t in df.Time]
    slope = zeros(length(df.Time))
    dt = 0.1
    t = df.Time[1]
    slope[1] = (predict(model, t+dt)-predict(model,t)) / dt
    t = df.Time[end]
    slope[end] = (predict(model, t)-predict(model,t-dt)) / dt
    for i=2:length(df.Time)-1
        t = df.Time[i]
        slope[i] = (predict(model, t+dt)-predict(model,t-dt)) / (2*dt)
    end
    return loess_height, slope
end

df = DataFrame(CSV.File("data/timelapses/database.csv"))
df.avg_height = abs.(df.avg_height)                      # Remove negative nums

forward_change = []
for st in unique(df.Strain)
    tf = filter(row->row.Strain.==st, df);
    for repli in unique(tf.Replicate)
        rtf = filter(row->row.Replicate.==repli, tf);
        Δh = d_height(rtf)
        append!(forward_change, Δh)
    end
end

df.forward_change = forward_change


##
chole = filter(row->row.Strain.=="SN503"&&row.Time.<=48&&row.Replicate.=="B", df);
loc_slo = local_slope(chole, 8.0)

p1 = scatter(chole.Time, chole.forward_change, label="Discrete")
p1 = scatter!(chole.Time, loc_slo, label="Local slope")

p2 = scatter(chole.avg_height, chole.forward_change, label="Discrete")
p2 = scatter!(chole.avg_height, loc_slo, label="Local slope")
plot(p1, p2, layout=(2,1))

##
plot(chole.Time, chole.avg_height)

##
##
time_window = 5.0
df = filter(row->row.Strain.=="SN503"&&row.Time.<=48&&row.Replicate.=="B", df);
per_window = time_window / (df.Time[end]-df.Time[1])
model = loess(df.Time, df.avg_height, span=per_window, degree=1)
loess_height = [predict(model, t) for t in df.Time]

scatter(reduce(vcat,my_keys)[2:end])
scatter!(df.Time, xlim=(0, 20), ylim=(0, 10))
##

##



h, s = smooth_heights(df, 4.0)
p1 = scatter(df.Time, [df.avg_height, h])
p2 = scatter(df.Time, [df.forward_change, s])
p3 = scatter(df.avg_height, df.forward_change)
p3 = scatter!(h, s)