using Plots, StatsPlots, Colors
using DataFrames, CSV

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

ecoli = filter(row->row.Strain.=="JT305L" &&row.Time.<=48, df);
aerom = filter(row->row.Strain.=="BGT127"&&row.Time.<=48, df);
pyeas = filter(row->row.Strain.=="pyeast"&&row.Time.<=48, df);
chole = filter(row->row.Strain.=="SN503"&&row.Time.<=48, df);

# Color definitions
ecoli_colors = colormap("Blues", 6)[3:5]'
aerom_colors = colormap("Oranges", 6)[3:5]'
pyeas_colors = colormap("Greens", 6)[3:5]'

## Panel A: sample biofilm slices

## Panel B: height dynamics for ecoli, aeromonas, pyeast
panel_b = @df ecoli scatter(:Time, :avg_height, group=:Replicate, label=false, alpha=0.9,color=ecoli_colors)
@df aerom scatter!(:Time, :avg_height, group=:Replicate, label=false,alpha=0.9, color=aerom_colors)
@df pyeas scatter!(:Time, :avg_height, group=:Replicate, label=false, alpha=0.9,color=pyeas_colors)
plot!(grid=false,xlabel="Time [h]", ylabel="Height [μm]", background_color=:transparent)

## Panel C: z' vs t
panel_c = @df ecoli scatter(:Time, :forward_change, group=:Replicate, label=false, alpha=0.9,color=ecoli_colors)
@df aerom scatter!(:Time, :forward_change, group=:Replicate, label=false,alpha=0.9, color=aerom_colors)
@df pyeas scatter!(:Time, :forward_change, group=:Replicate, label=false, alpha=0.9,color=pyeas_colors)
plot!(grid=false,xlabel="Time [h]", ylabel="ΔHeight [μm/h]", background_color=:transparent)

## Panel D: z' vs z 
panel_d = @df ecoli scatter(:avg_height, :forward_change, group=:Replicate, label=false, alpha=0.9,color=ecoli_colors)
@df aerom scatter!(:avg_height, :forward_change, group=:Replicate, label=false,alpha=0.9, color=aerom_colors)
@df pyeas scatter!(:avg_height, :forward_change, group=:Replicate, label=false, alpha=0.9,color=pyeas_colors)
plot!(grid=false,xlabel="Height [μm]", ylabel="ΔHeight [μm/h]", background_color=:transparent)

## All together
lay = @layout [a{0.5w} [b
                        c]]
plot(panel_b, panel_c, panel_d, layout=lay)



## Exploring height changes.
# Loop over individual timelapses
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
function d_height(df)
    h_change = zeros(size(df)[1]) 
    h_change[end] = NaN
    if size(df)[1] > 1
        h_change[1:end-1] = (df.avg_height[2:end]-df.avg_height[1:end-1]) ./ 
                            (df.Time[2:end]-df.Time[1:end-1])
    end
    return h_change
end

d_height(pyeas)

##
@df chole scatter(:Time, :forward_change, group=:Replicate, color=[1 2 3])
##
@df df scatter(:Time, :mid_height, group=(:Replicate, :Strain), xlim=(0, 48), alpha=0.5, legend=:topleft, legendfontsize=4)
