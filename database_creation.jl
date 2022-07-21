#=
This script loops over the raw profiles + cleaning data on /data/timelapses
calculates relevant scalar metrics and then saves all into a single 
file for plotting. This database.csv output contains only experimental 
outputs, no models included. 
=#
using DataFrames, NPZ, CSV
using Statistics, NaNMath, LsqFit, Glob

function d_height(df)
    h_change = zeros(size(df)[1]) 
    h_change .= NaN
    if size(df)[1] > 2
        h_change[1:end-1] = (df.avg_height[2:end]-df.avg_height[1:end-1]) ./ 
                            (df.time[2:end]-df.time[1:end-1])
    end
    return h_change
end

function smooth_heights(df, dt)
    model(x, p) = p[1] .+ p[2]*x # Linear model
    y_smooth = zeros(size(df)[1])
    slope_mean = zeros(size(df)[1])
    slope_error = zeros(size(df)[1])
    y_smooth .= NaN              # Start them as NaNs and then fill
    slope_mean .= NaN
    slope_error .= NaN
    x, y = df.time, df.avg_height
    if size(df)[1] > 2           # If we have more points!
        for i=1:length(x)
            idx = (x .> x[i]-dt/2) .&& (x .< x[i]+dt/2)
            x_c, y_c = x[idx], y[idx]
            p_guess = [y_c[1], (y_c[end]-y_c[1])/(x_c[end]-x_c[1])]
            fit = curve_fit(model, x_c, y_c, p_guess)
            y_smooth[i] = model(x[i], fit.param)
            slope_mean[i] = fit.param[2]
            slope_error[i] = sqrt(estimate_covar(fit)[2,2])
        end
    end
    return y_smooth, slope_mean, slope_error
end

"""
Takes a folder as input, and outputs a 
DataFrame containing profiles and some
simple metrics
"""
function add_to_database(folder)
    strain = match(r"(\d{4})-(\d{2})-(\d{2})_(\w+)", folder)[4]
    timelapses = glob("profiles_*.npy", folder*"/")
    df = DataFrame(CSV.File(folder*"/cleaning.csv"))
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
    for n in findall(isnan.(df.avg_height))     # Fix NaN values with linear interpolation, print a warning.
        print("NaN height found in file:"*df.file[n])
        df.avg_height[n] = 0.5*df.avg_height[n-1]+0.5*df.avg_height[n+1]
        df.std_height[n] = 0.5*df.std_height[n-1]+0.5*df.std_height[n+1]
    end
    forward_change = []
    smooth_height = []
    slope = []
    slope_error = []
    order = []
    for repli in unique(df.replicate)
        tf = filter(row->row.replicate.==repli, df);
        dh = d_height(tf)
        h, s, se = smooth_heights(tf, 4.0)
        my_order = Array(1:size(tf)[1])
        append!(forward_change, dh)
        append!(smooth_height, h)
        append!(slope, s)    
        append!(slope_error, se)  
        append!(order, my_order)    
    end
    df.forward_change = forward_change 
    df.smooth_height = smooth_height
    df.slope = slope
    df.slope_error = slope_error
    df.order = order
    return df[!, Not(:profile)]
end

rt_folder = "data/timelapses"
tl_folders = readdir(rt_folder)[1:10]
Data = DataFrame(file = String[], replicate = String[], time=Float32[], border_l = Int[],
                 border_r = Int[], offset_l = Float32[], offset_r = Float32[], 
                 displacement = Int32[], zoom =Int32[], strain = String[], mid_height = Float32[],
                 max_height = Float32[], min_height = Float32[], width = Float32[], hL = Int[], hR = Int[],
                 avg_height = Float32[], std_height = Float32[], forward_change = Float32[], 
                 smooth_height = Float32[], slope = Float32[], slope_error= Float32[], order=Int32[])
##
for folder in tl_folders
    println(folder)
    strain_dataframe = add_to_database(rt_folder*"/"*folder*"/")
    strain_dataframe = strain_dataframe[(strain_dataframe.border_l .!= 0), :]
    calculated = calculations(strain_dataframe)
    global Data = vcat(Data, calculated)
end
##
# Round some numbers up
Data.time = round.(Data.time, digits=2)
Data.offset_l = round.(Data.offset_l, digits=1)
Data.offset_r = round.(Data.offset_r, digits=1)
Data.mid_height = round.(Data.mid_height, digits=3)
Data.max_height = round.(Data.max_height, digits=3)
Data.min_height = round.(Data.min_height, digits=3)
Data.width = round.(Data.width, digits=3)
Data.avg_height = round.(Data.avg_height, digits=3)
Data.std_height = round.(Data.std_height, digits=3)
Data.forward_change = round.(Data.forward_change, digits=3)
Data.smooth_height = round.(Data.smooth_height, digits=3)
Data.slope = round.(Data.slope, digits=3)
Data.slope_error = round.(Data.slope_error, digits=3)
##
CSV.write("data/timelapses/database.csv", Data)
##
