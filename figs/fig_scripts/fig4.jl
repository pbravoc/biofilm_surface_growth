using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Colors
using LsqFit

function get_average(df, strain_name)
    tf = filter(x->x.strain .== strain_name,df)
    l = Int(size(tf)[1]/3)
    h = reshape(tf.avg_height, (l, 3))
    h_avg = reduce(vcat, mean(h, dims=2))
    h_std = reduce(vcat, std(h, dims=2))
    t = tf.time[l+1:2*l]
    scatter!(t, h_avg, ribbon=h_std, alpha=0.3,
             fillalpha=0.6, label=strain_name,
             markersize=3, linewidth=2)
end

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
df = filter(x->x.replicate in ["A", "B", "C"] && x.time .< 48, Df)
df2 = filter(x->x.replicate in unique(Df.replicate)[4:end], Df)
@df df scatter(:time, :avg_height, markersize=1, legend=false)
allfit = DataFrame(CSV.File("data/sims/allpointsol.csv"))
c0, c1, c2, c3 = :black, 1, 2, 3
#c0, c1, c2, c3 = :black, :black, :black, :black

# Heights vs time
p1_1 = @df filter(x->x.strain == "bgt127", df2) scatter(:time, :avg_height, yerror=:std_height,color=c0, markersize=3,  label="Aeromonas")
bgt127_lt = DataFrame(CSV.File("data/sims/bgt127_lt.csv"))
@df bgt127_lt plot!(:time, :interface, color=c1, linestyle=[:dash :dot], linewidth=2, label=false)
@df allfit plot!(:t, :bgt127, color=c1, linewidth=2, label=false)
p1_2 = @df filter(x->x.strain == "jt305", df2) scatter(:time, :avg_height, yerror=:std_height, color=c0, markersize=3, label="E coli")
jt305_lt = DataFrame(CSV.File("data/sims/jt305_lt.csv"))
@df jt305_lt plot!(:time, :interface, color=c2, linestyle=[:dash :dot], linewidth=2, label=false)
@df allfit plot!(:t, :jt305, color=c2, linewidth=2, label=false)
p1_3 = @df filter(x->x.strain == "gob33", df2) scatter(:time, :avg_height, yerror=:std_height,color=c0,  markersize=3, label="Yeast (aa)")
plot!(legend=:topleft, xlabel="Time [hr]", ylabel="Height [μm]")
gob33_lt = DataFrame(CSV.File("data/sims/gob33_lt.csv"))
@df gob33_lt plot!(:time, :interface, color=c3, linestyle=[:dash :dot], linewidth=2, label=false)
@df allfit plot!(:t, :gob33, color=c3, linewidth=2, label=false)
p1 = plot(p1_1, p1_2, p1_3, layout=(1,3), size=(900, 250), legend=false, ylim=(0, 800), xlabel="Time [hr]", ylabel="Height [μm]")
plot!(bottom_margin=4mm, left_margin=4mm)
# Height distributions
pf = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
p2 = plot()
pf2 = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
pf2.h_max = pf2.x1 .*pf2.x3 ./ pf2.x2
p2_1 = density(sort(pf[pf.strain.=="bgt127",4])[6:974], bins=20, xticks=[206, 240], color=c1, fill=:true, fillalpha=0.3)
temp_h = filter(x->x.fit in ["all", "48h"] && x.strain .== "bgt127", pf2)
scatter!(temp_h[:,6], [0.0835, 0.0835], color=c1, marker=[:circle, :diamond], legend=false)
annotate!(217, 0.083, "223.8μm", 6)
annotate!(231, 0.083, "224.5μm", 6)
plot!(yticks=false)
p2_2 = density(sort(pf[pf.strain.=="jt305",4])[6:974], color=c2, xticks=[90, 520],fill=:true, fillalpha=0.3)
temp_h = filter(x->x.fit in ["all", "48h"] && x.strain .== "jt305", pf2)
scatter!(temp_h[:,6], [0.0016, 0.0067], color=c2, marker=[:circle, :diamond], legend=false)
annotate!(410, 0.0016, "307.1μm", 6)
annotate!(315, 0.0067, "213.6μm", 6)
plot!(yticks=false)
p2_3 = density(sort(pf[pf.strain.=="gob33",4])[6:974], color=c3, fill=:true, fillalpha=0.3, xticks=[295, 1700])
temp_h = filter(x->x.fit in ["all", "48h"] && x.strain .== "gob33", pf2)
scatter!(temp_h[:,6], [0.0011, 0.003], color=c3, marker=[:circle, :diamond], legend=false)
annotate!(1100, 0.0011, "696.3μm", 6)
annotate!(900, 0.003, "471.6μm", 6)
plot!(yticks=false)
p2 = plot(p2_1, p2_2, p2_3, layout=(1,3), legend=false, size=(700, 250))
# Residuals
df = DataFrame(CSV.File("data/timelapses/fit_predictions.csv"))
df[7,4] = 0.4
p3_1 = @df df[df.strain .== "bgt127",:] plot(:t, :h-:h, ribbon=:h_std,linewidth=2, alpha=0, color=c0, fillalpha=0.15, label=false)
@df df[df.strain .== "bgt127",:] plot!(:t, [:interface_48-:h, :interface_all-:h], linewidth=2, color=c1, linestyle=[:solid :dash], label=false)
p3_2 = @df df[df.strain .== "jt305",:] plot(:t, :h-:h, ribbon=:h_std,linewidth=2, alpha=0, color=c0, fillalpha=0.15, label=false)
@df df[df.strain .== "jt305",:] plot!(:t, [:interface_48-:h, :interface_all-:h], linewidth=2, color=c2, linestyle=[:solid :dash], label=false)
p3_3 = @df df[df.strain .== "gob33",:] plot(:t, :h-:h, ribbon=:h_std,linewidth=2, alpha=0, color=c0, fillalpha=0.15, label=false)
@df df[df.strain .== "gob33",:] plot!(:t, [:interface_48-:h, :interface_all-:h], linewidth=2, color=c3, linestyle=[:solid :dash], label=false)

l = @layout [ grid(1,3){0.8h}
              b{0.2333w} c{0.1w} d{0.2333w} e{0.1w} f{0.2333w} g{0.1w}]
plot(p1_1, p1_2, p1_3, 
     p3_1, p2_1, 
     p3_2, p2_2, 
     p3_3, p2_3, layout=l, size=(900, 350))
##

##
plot(p1, p2, p3, p4, size=(900, 450), layout=l, grid=false)
#savefig("figs/fig4/fig4_notfinal.pdf")
##