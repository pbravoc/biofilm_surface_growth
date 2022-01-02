using Plots, StatsPlots, Colors
using DataFrames, CSV

df = DataFrame(CSV.File("data/timelapses/database.csv"))    # Load dataset
#df.avg_height = abs.(df.avg_height)                      # Remove negative nums
#tf = filter(row->row.Strain∈["JT305L", "BGT127","pyeast"]&& # Still ugly, this goes in
tf = filter(row->row.Strain∈["SN503", "JT1080","BGT127"]&& # Pretty!
                    row.Replicate∈["A","B","C"]&&
                    row.Time.<=48, df);                 # Select only 3 strains
my_colors = [colormap("Blues", 6)[3:5]                  # Color definitions
             colormap("Oranges", 6)[3:5]
             colormap("Greens", 6)[3:5]]'

## Panel A: sample biofilm slices

panel_b = @df tf scatter(:Time, :loess_height, group=(:Strain, :Replicate), 
                         legend=:topleft, alpha=0.5,color=my_colors)
plot!(grid=false,xlabel="Time [h]", ylabel="Height [μm]", 
      background_color=:transparent)
panel_c = @df tf scatter(:Time, :local_slope, group=(:Strain, :Replicate), 
                         label=false, alpha=0.5,color=my_colors)
plot!(grid=false,xlabel="Time [h]", ylabel="ΔHeight [μm/h]", 
      background_color=:transparent)
panel_d = @df tf scatter(:loess_height, :local_slope, group=(:Strain, :Replicate), 
                         label=false, alpha=0.5,color=my_colors)
plot!(grid=false,xlabel="Height [μm]", ylabel="ΔHeight [μm/h]", 
      background_color=:transparent)
lay = @layout [a{0.5w} [b
                        c]]
plot(panel_b, panel_c, panel_d, layout=lay)     ## All together

