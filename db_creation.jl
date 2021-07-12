using DataFrames, NPZ 
using Arrow

"""
Takes a dictionary as input, and outputs a 
DataFrame containing profiles and some
simple metrics
"""
function add_to_database(df, dict)
    times = npzread(string(dict["folder"], "times.npy"))
    replicate = ["A", "B", "C"]
    for i=1:3
        profile = npzread(string(dict["folder"], "profiles_",    # Open current working replicate
                          replicate[i],".npy"))
        for j=1:size(times)[2]
            rowdata = (dict["strain"], dict["date"], replicate[i], times[j], j,
                       dict["zoom"], dict["by"], profile[j,:])
            push!(df, rowdata)
        end
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