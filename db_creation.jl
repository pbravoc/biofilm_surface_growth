using DataFrames, NPZ 
using Arrow

"""
Takes a dictionary as input, and outputs a 
DataFrame containing profiles and some
simple metrics
"""
function add_to_database(df, dict)
    times = npzread(string(dict["folder"], "times.npy"))
    replicate = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]

    # Main (timelapse) points
    for i=1:3
        profile = npzread(string(dict["folder"], "profiles_",    
                          replicate[i],".npy"))
        for j=1:size(times)[2]
            rowdata = (dict["strain"], dict["date"], replicate[i], times[i,j], j,
                       dict["zoom"], dict["by"], profile[j,:])
            push!(df, rowdata)
        end
    end

    # Secondary (control) points 
    control_times = npzread(string(dict["folder"], "times_control.npy"))
    control_profiles = npzread(string(dict["folder"], "profiles_control.npy"))
    for i=1:size(control_times)[1]
        zoom = i < 4 ? 50 : 10
        rowdata = (dict["strain"], dict["date"], replicate[i+3], control_times[i], i,
        zoom, dict["by"], control_profiles[i,:])
        push!(df, rowdata)
    end
end

# Initialize an empty database
df = DataFrame(Strain = String[], Date = String[], Replicate = String[], 
               Time = Float32[], Order=Int32[], Zoom = Float32[], 
               By = String[], Profile = Array[]);

# BGT127: Aeromonas
bgt127 = Dict("folder" => "data/timelapses/2021-06-25_bgt127/",
              "strain" => "BGT127", "date" => "2021-06-25", 
              "zoom" => 50, "by" => "pbravo")

# JT305: Ecoli
jt305 = Dict("folder" => "data/timelapses/2021-07-09_jt305/",
              "strain" => "JT305", "date" => "2021-07-09", 
              "zoom" => 50, "by" => "pbravo")

add_to_database(df, bgt127)
#add_to_database(df, jt305)


# Middle_height
df.mid_height = [df.Profile[i][Int(length(df.Profile[i])/2)] for i=1:size(df)[1]]



Arrow.write("data/timelapses/profile_database.arrow", 
            df , compress = :zstd)