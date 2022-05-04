#= This code loops over the database, and returns the best fits
for the 48 hour data. One for each timelapse, and one for the aggregated
data.
It also calculated the best fit using 48h data + measurements from 
longtime_data.csv
=#
include("surface_functions.jl")

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
Df = filter(x-> x.avg_height .> 0, Df)                    # Smaller than 0 values don't make physical sense
df2 = DataFrame(CSV.File("data/timelapses/longtime_data.csv"))
strain_list = unique(Df.strain)
long_list = ["bgt127", "gob33", "jt305"]
P = []
Strain = []
Fit = []
bounded = []
## Details on each model to fit
logistic_n = model_struct([0.1, 1.0],
                          [0.8, 0.1, 0.9, 0.12],
                          [0.1, 0.01, 0.01, 0.001],
                          [5.0, 1e1, 1.0, 1.0],
                           ode_logistic_n, "logistic_n")
##
for model_choice in [logistic_n]
    # Get the best fits for less than 48h 
    for bound in  [false, true]
        for strain in strain_list 
            df = filter(x-> x.strain .== strain && x.time .<48, Df)              
            fit_params = fit_data(df.time, df.avg_height, model_choice, bound)
            append!(P, [fit_params.u])
            append!(Strain, [strain])
            append!(Fit, ["48h"])
        end
        # Get the best fits for each timelapse
        for strain in strain_list 
            for replicate in ["A", "B", "C"]     
                df = filter(x-> x.strain .== strain && x.time .<48 &&
                                x.replicate .== replicate, Df)
                fit_params = fit_data(df.time, df.avg_height, model_choice, bound)
                append!(P, [fit_params.u])
                append!(Strain, [strain])
                append!(Fit, [replicate])
            end
        end
        # This is to get the 'long time' best fit.
        if strain in long_list
            df = filter(x-> x.strain in strain_list && x.time .< 48 &&
                            x.replicate in ["A", "B", "C"], Df)
            tf = filter(x-> x.strain .== strain, df)
            tf2 = filter(x-> x.strain .== strain, df2)
            t = append!(tf.time, tf2.time)
            h = append!(tf.avg_height, tf2.avg_height)
            fit_params = fit_data(t, h, model_choice, bound)
            append!(P, [fit_params.u])
            append!(Strain, [strain])
            append!(Fit, ["long"])
        end
        append!(bounded, repeat([bound], length(Strain)))
    end
    #
    pf = hcat(DataFrame("strain"=>Strain, "fit"=>Fit, "bounded"=>bounded),
            DataFrame(Matrix(reduce(hcat, P)'), :auto))
    ## Save to file
    CSV.write("data/timelapses/fit_params_"*model_choice.name*".csv", pf)
end
##
model_choice = logistic_n
for strain in strain_list 
    df = filter(x-> x.strain .== strain && x.time .<48, Df)       
    print(length(df.time))       
    fit_params = fit_data(df.time, df.avg_height, model_choice, true)
    append!(P, [fit_params.u])
    append!(Strain, [strain])
    append!(Fit, ["48h"])
end