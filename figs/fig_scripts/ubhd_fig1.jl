using Plots, StatsPlots, Colors
using DataFrames, CSV

df = DataFrame(CSV.File("data/timelapses/database.csv"))
df.avg_height = abs.(df.avg_height)                      # Remove negative nums
ecoli = filter(row->row.Strain.=="JT305L" &&row.Time.<=48, df);
aerom = filter(row->row.Strain.=="BGT127"&&row.Time.<=48, df);
pyeas = filter(row->row.Strain.=="pyeast"&&row.Time.<=48, df);

## Color definitions

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

panel_c = @df ecoli scatter(:avg_height, group=:Replicate, label=false, alpha=0.9,color=ecoli_colors)
@df aerom scatter!(:avg_height, group=:Replicate, label=false,alpha=0.9, color=aerom_colors)
@df pyeas scatter!(:avg_height, group=:Replicate, label=false, alpha=0.9,color=pyeas_colors)
plot!(grid=false,xlabel="Time [h]", ylabel="Height [μm]", background_color=:transparent)

## Panel D: z' vs z 
@df pyeas scatter(:Time, :avg_height, group=:Replicate, legend=:topleft, color=[1 2 3])

## All together
