using DataFrames, NPZ, Arrow, CSV
using Statistics, NaNMath, Loess, Glob

function d_height(df)
    h_change = zeros(size(df)[1]) 
    h_change .= NaN
    if size(df)[1] > 2
        h_change[1:end-1] = (df.avg_height[2:end]-df.avg_height[1:end-1]) ./ 
                            (df.time[2:end]-df.time[1:end-1])
    end
    return h_change
end

function smooth_heights(df, time_window)
    loess_height, slope = zeros(size(df)[1]), zeros(size(df)[1])
    loess_height .= NaN
    slope .= NaN
    if size(df)[1] > 2
        per_window = time_window / (df.time[end]-df.time[1])
        model = loess(df.time, df.avg_height, span=per_window, degree=1)
        loess_height = [predict(model, Float64(df.time[i])) for i=1:length(df.time)]
        slope = zeros(length(df.time))
        dt = 0.1
        t = df.time[1]
        slope[1] = (predict(model, Float64(t+dt))-predict(model,Float64(t))) / dt
        t = df.time[end]
        slope[end] = (predict(model, Float64(t))-predict(model,Float64(t-dt))) / dt
        for i=2:length(df.time)-1
            t = df.time[i]
            slope[i] = (predict(model, Float64(t+dt))-predict(model,Float64(t-dt))) / (2*dt)
        end
    end
    return loess_height, slope
end

"""
Takes a folder as input, and outputs a 
DataFrame containing profiles and some
simple metrics
"""
function add_to_database(folder)
    strain = match(r"(\d{4})-(\d{2})-(\d{2})_(\w+)", folder)[4]
    timelapses = glob("profiles_*.npy", folder*"/Clean/")
    df = DataFrame(CSV.File(folder*"/Clean/cleaning.csv"))
    profile = []
    for tl in timelapses
        curr_profile = npzread(tl)
        for i in range(1, size(curr_profile)[1])
            append!(profile, [curr_profile[i,:]])
        end
    end
    df.profile = profile
    df.strain = repeat([strain], size(df)[1])
    return df
end

function calculations(df)
    df.mid_height = [df.profile[i][Int(floor(length(df.profile[i])/2))] for i=1:size(df)[1]]
    df.max_height = [NaNMath.maximum(df.profile[i]) for i=1:size(df)[1]]
    df.min_height = [NaNMath.minimum(df.profile[i]) for i=1:size(df)[1]]
    l = [findall(x->!isnan(x), y)[1] for y in df.profile]
    r = [findall(x->!isnan(x), y)[end] for y in df.profile]
    df.width = (r-l)* 0.17362 * 1e-3 * 50 ./ df.zoom
    h = 2e3 ./ (0.17362 * 50 ./ df.zoom)
    df.hL = Int.(floor.(length.(df.profile)/2 .- h./2))
    df.hR = Int.(floor.(length.(df.profile)/2 .+ h./2))
    df.avg_height = [NaNMath.mean(df.profile[i][df.hL[i]:df.hR[i]]) for i=1:size(df)[1]]
    df.std_height = [NaNMath.std(df.profile[i][df.hL[i]:df.hR[i]]) for i=1:size(df)[1]]
    forward_change = []
    loess_height = []
    local_slope = []
    for repli in unique(df.replicate)
        tf = filter(row->row.replicate.==repli, df);
        dh = d_height(tf)
        h, s = smooth_heights(tf, 4.0)
        append!(forward_change, dh)
        append!(loess_height, h)
        append!(local_slope, s)    
    end
    df.forward_change = forward_change 
    df.loess_height = loess_height
    df.local_slope = local_slope
    return df[!, Not(:profile)]
end

##
rt_folder = "/run/media/pablo/T7/Documents/Research/Biofilms/Data/Interferometry/radial_timelapses/"
tl_folders = readdir(rt_folder)[1:6]
Data = DataFrame()
for folder in tl_folders
    strain_dataframe = add_to_database(rt_folder*folder*"/")
    calculated = calculations(strain_dataframe)
    #append!(Data, calculated)
    global Data = vcat(Data, calculated)
end
##
CSV.write("data/timelapses/database_updated.csv", Data)

##

