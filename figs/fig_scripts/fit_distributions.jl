using DataFrames, CSV
using Plots, StatsPlots

df = DataFrame(CSV.File("data/sims/bootstrap/allfits.csv"))
##
@df df plot(:α, :β, :L, grouping=:strain)
##
@df df histogram(:α, grouping=:strain)
##
function plot_variable(df, strain_name, c)
    tf = filter(x-> x.strain .== strain_name, df)
    p1 = @df tf density(:α, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="α")
    vline!([tf[tf.name .== "fit_48",:].α, tf[tf.name .== "fit_all",:].α], 
           color=:black, linestyle=[:dash :solid], legend=false)
    p2 = @df tf density(:β, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="β")
    vline!([tf[tf.name .== "fit_48",:].β, tf[tf.name .== "fit_all",:].β], 
            color=:black, linestyle=[:dash :solid], legend=false) 
    p3 = @df tf density(:L, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="L")
    vline!([tf[tf.name .== "fit_48",:].L, tf[tf.name .== "fit_all",:].L], 
            color=:black, linestyle=[:dash :solid], legend=false)
    p4 = @df tf density(:h_max, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="h_max")
    vline!([tf[tf.name .== "fit_48",:].h_max, tf[tf.name .== "fit_all",:].h_max], 
            color=:black, linestyle=[:dash :solid], legend=false)     
    plot(p1, p2, p3, p4, layout=(4,1))
    return [p1, p2, p3, p4]
end
##
p = plot_variable(df, "gob33", 2)
p1 = plot(p[1])
p2 = plot(p[2])
p3 = plot(p[3])
p4 = plot(p[4], xlim=(0, 2500))
plot(p1, p2, p3, p4)
##
function plot_variable_outlier(df, strain_name, c)
        tf = filter(x-> x.strain .== strain_name, df)
        α = sort(tf.α)[26:974]
        β = sort(tf.β)[26:974]
        L = sort(tf.L)[26:974]
        h_max = sort(tf.h_max)[26:974]
        print(h_max[end])
        p1 = @df tf density(α, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="α [μm/hr]")
        vline!([tf[tf.name .== "fit_48",:].α, tf[tf.name .== "fit_all",:].α], 
               color=:black, linestyle=[:dash :solid], legend=false)
        p2 = @df tf density(β, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="β [μm/hr]")
        vline!([tf[tf.name .== "fit_48",:].β, tf[tf.name .== "fit_all",:].β], 
                color=:black, linestyle=[:dash :solid], legend=false) 
        p3 = @df tf density(L, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="L [μm]")
        vline!([tf[tf.name .== "fit_48",:].L, tf[tf.name .== "fit_all",:].L], 
                color=:black, linestyle=[:dash :solid], legend=false)
        p4 = @df tf density(h_max, fillrange = 0, color=c, fillalpha = 0.5, fillcolor = c, title="h_max [μm]")
        vline!([tf[tf.name .== "fit_48",:].h_max, tf[tf.name .== "fit_all",:].h_max], 
                color=:black, linestyle=[:dash :solid], legend=false)     
        plot(p1, p2, p3, p4, layout=(4,1))
        return [p1, p2, p3, p4]
end
p = plot_variable(df, "bgt127", 1)
p1 = plot(p[1])
p2 = plot(p[2])
p3 = plot(p[3])
p4 = plot(p[4])
plot(p1, p2, p3, p4, size=(800, 500))
savefig("figs/bootstrapping/bgt127.svg")

##
function remove_outlier(df, strain_name)
        tf = filter(x-> x.strain .== strain_name, df)

        idxs = sortperm(tf.h_max)[26:974]
        return tf.α[idxs], tf.β[idxs], tf.L[idxs] 
end
bgt127_data = remove_outlier(df, "bgt127")
jt305_data = remove_outlier(df, "jt305")
gob33_data = remove_outlier(df, "gob33")
##
scatter(bgt127_data[1], bgt127_data[2], bgt127_data[3], alpha=1.0,markersize=1, markerstrokealpha=0.0, label="BGT127-Aeromonas")
scatter!(jt305_data[1], jt305_data[2], jt305_data[3], alpha=1.0,markersize=1, markerstrokealpha=0.0, label="JT305-Ecoli")
scatter!(gob33_data[1], gob33_data[2], gob33_data[3], alpha=1.0,markersize=1, markerstrokealpha=0.0, label="GOB33-Petite Yeast")
plot!(xlabel="Alpha [um/hr]", ylabel="Beta [um/hr]", zlabel="L [um]", size=(1000,1000))
savefig("figs/bootstrapping/scatter.html")
##
scatter(bgt127_data[1], bgt127_data[2], alpha=0.8, label="BGT127")
scatter!(jt305_data[1], jt305_data[2], alpha=0.8, label="JT305")
scatter!(gob33_data[1], gob33_data[2], alpha=0.8, label="GOB33")
plot!(xlabel="α [μm/hr]", ylabel="β [μm/hr]", legend=:bottomright)
##
scatter(bgt127_data[1], bgt127_data[3], alpha=0.8, label="BGT127")
scatter!(jt305_data[1], jt305_data[3], alpha=0.8, label="JT305")
scatter!(gob33_data[1], gob33_data[3], alpha=0.8, label="GOB33")
plot!(xlabel="α [μm/hr]", ylabel="L [μm]")

##
scatter(bgt127_data[2], bgt127_data[3], alpha=0.8, label="BGT127")
scatter!(jt305_data[2], jt305_data[3], alpha=0.8, label="JT305")
scatter!(gob33_data[2], gob33_data[3], alpha=0.8, label="GOB33")
plot!(xlabel="β [μm/hr]", ylabel="L [μm]")

##
marginalkde(bgt127_data[1], bgt127_data[2])
marginalkde!(jt305_data[1], jt305_data[2])

#1. Set the goal on slide 2
#2. We're not interested in bulk and average over the small distance.
#3. Exponential growth + decreased (go back to the left figure)
#4. Logistic model is bad (!= doesnt do so well)
#5. Steady state vs equilibrium