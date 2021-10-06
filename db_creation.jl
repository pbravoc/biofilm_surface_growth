using DataFrames, NPZ, Arrow, JLD2
using Statistics, NaNMath

"""
Takes a dictionary as input, and outputs a 
DataFrame containing profiles and some
simple metrics
"""
function add_to_database(df, dict, n, long, cont)
    times = npzread(string(dict["folder"], "times.npy"))
    replicate = ["A", "B", "C", "D", "E", "F", "G", "H", "I"] # Improve this
    repl_long = ["A", "B", "C", "D", "D", "E", "E", "F", "F", "G", "G", "H", "H", "I", "I", "J", "K", "L"]
    # Main (timelapse) points
    for i=1:n
        profile = npzread(string(dict["folder"], "profiles_",    
                          replicate[i],".npy"))
        if long 
            for j=1:(size(times)[2]-1)
                rowdata = (dict["strain"], dict["date"], replicate[i], times[i,j], j,
                        dict["zoom"], dict["by"], profile[j,:])
                push!(df, rowdata)
            end 
        else 
            for j=1:(size(times)[2])
                rowdata = (dict["strain"], dict["date"], replicate[i], times[i,j], j,
                        dict["zoom"], dict["by"], profile[j,:])
                push!(df, rowdata)
            end 
        end
    end

    if cont
        # Secondary (control) points 
        control_times = npzread(string(dict["folder"], "times_control.npy"))
        control_profiles = npzread(string(dict["folder"], "profiles_control.npy"))

        if long
            for i=1:size(control_times)[1]
                if i in [1,2,3,4,6,8]
                    zoom = 50
                else 
                    zoom = 10
                end
                rowdata = (dict["strain"], dict["date"], repl_long[i], control_times[i], i, zoom, dict["by"], control_profiles[i,:])
                push!(df, rowdata)
            end
        else 
            for i=1:size(control_times)[1]
                zoom = i < 4 ? 50 : 10           # ABCDEF are 50x, GHI are 10x due to size
                rowdata = (dict["strain"], dict["date"], replicate[i+3], control_times[i], i, zoom, dict["by"], control_profiles[i,:])
                push!(df, rowdata)
            end
        end
    end

    println(string("Added " , dict["strain"]))
end

# Initialize an empty database
df = DataFrame(Strain = String[], Date = String[], Replicate = String[], 
               Time = Float32[], Order=Int32[], Zoom = Float32[], 
               By = String[], Profile = Array[]);
println("Columns created")

# Add Vibrios
#Df = DataFrame(Arrow.Table("/home/pablo/Biofilms/Data/radialv2.arrow"));   # Lab pc
Df = DataFrame(Arrow.Table("/home/pablo/Documents/radialv2.arrow"));        # Omilen
tf = select(Df, :Strain, :Date, :Replicate, :Time, :Order, :Zoom, :By, :Profile);
for i=1:size(tf)[1]
    push!(df, tf[i,:])
end
println("Added older vcholerae data")

# Dictionaries for each strain 
# BGT127: Aeromonas
bgt127 = Dict("folder" => "data/timelapses/2021-06-25_bgt127/",
              "strain" => "BGT127", "date" => "2021-06-25", 
              "zoom" => 50, "by" => "pbravo")

# JT305: Ecoli
jt305 = Dict("folder" => "data/timelapses/2021-07-09_jt305/",
              "strain" => "JT305", "date" => "2021-07-09", 
              "zoom" => 50, "by" => "pbravo")

# JT305-long: Ecoli
jt305l = Dict("folder" => "data/timelapses/2021-08-27_jt305/",
              "strain" => "JT305L", "date" => "2021-08-27", 
              "zoom" => 50, "by" => "pbravo")

# Yeast (LB)
yeast = Dict("folder" => "data/timelapses/2021-07-23_yeast/",
              "strain" => "yeast", "date" => "2021-07-23", 
              "zoom" => 50, "by" => "pbravo")

# Bacillus
bacillus = Dict("folder" => "data/timelapses/2021-07-30_bacillus/",
              "strain" => "bacillus", "date" => "2021-07-30", 
              "zoom" => 50, "by" => "pbravo")

# Petite yeast
pyeast = Dict("folder" => "data/timelapses/2021-09-03_pyeast/",
              "strain" => "pyeast", "date" => "2021-09-03", 
              "zoom" => 50, "by" => "pbravo")

# Add metadata and profiles to database 
add_to_database(df, bgt127, 3, false, true)
add_to_database(df, jt305, 3, false, true)
add_to_database(df, jt305l, 3, true, true)
#add_to_database(df, yeast, 3, false, false)
#add_to_database(df, bacillus, 2, false, false)
#add_to_database(df, pyeast, 3, false, false)

# Change infinites for NaNs 
for x in df.Profile
    x[isinf.(x)] .= NaN 
end

# Simple calculations
df.mid_height = [df.Profile[i][Int(floor(length(df.Profile[i])/2))] for i=1:size(df)[1]]
df.max_height = [NaNMath.maximum(df.Profile[i]) for i=1:size(df)[1]]
l = [findall(x->!isnan(x), y)[1] for y in df.Profile]
r = [findall(x->!isnan(x), y)[end] for y in df.Profile]
df.width = (r-l)* 0.17362 * 1e-3 * 50 ./ df.Zoom

h = 2e3 ./ (0.17362 * 50 ./ df.Zoom)
df.hL = Int.(floor.(length.(df.Profile)/2 .- h./2))
df.hR = Int.(floor.(length.(df.Profile)/2 .+ h./2))
df.avg_height = [NaNMath.mean(df.Profile[i][df.hL[i]:df.hR[i]]) for i=1:size(df)[1]]
df.std_height = [NaNMath.std(df.Profile[i][df.hL[i]:df.hR[i]]) for i=1:size(df)[1]]
#df.volume  
#df.std

# Complex calculations 
#df.Roughness 
#df.Fractal 
#df.Curvature

# Remove profiles from database to make it small
df = select(df, Not(:Profile))

# Write dataframe as an arrow file
jldsave("data/timelapses/profile_database.jld2"; df)

#Arrow.write("data/timelapses/profile_database.arrow", 
#            df)#, compress = :zstd)
println("success!")