# This code is ugly, but since data is kind of different and i'm
# being lazy it works. Write a reasonable script later?

using DataFrames, CSV
using Statistics, NaNMath
using Plots, StatsPlots, ColorSchemes, Plots.Measures

Df =  DataFrame(CSV.File("data/timelapses/database.csv"))
final_heights = []
heights_error = []
comparison_time = []
# bgt127
tf = filter(x->x.strain .== "bgt127", Df)
final_index = [100, 200, 300, 301, 302, 303]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:6],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                         std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# jt305
tf = filter(x->x.strain .== "jt305", Df)
final_index = [128, 256, 384, 385, 387, 389]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:6],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# gob33
tf = filter(x->x.strain .== "gob33", Df)
final_index = [102, 204, 306, 307, 308, 309]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:6],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# y55
tf = filter(x->x.strain .== "y55", Df)
final_index = [144, 288, 432, 433, 434, 435]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:6],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# bh1514
tf = filter(x->x.strain .== "bh1514", Df)
final_index = [83, 166, 249, 250, 251]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:5],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:5],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# ea387
tf = filter(x->x.strain .== "ea387", Df)
final_index = [82, 164, 246, 247 , 248, 249]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:5],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))

# cc117
tf = filter(x->x.strain .== "cc117", Df)
final_index = [88, 176, 264, 265 , 266, 267]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:5],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# sw520
tf = filter(x->x.strain .== "sw520", Df)
final_index = [78, 156, 234, 235 , 236, 237]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:5],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
# sw519
tf = filter(x->x.strain .== "sw519", Df)
final_index = [89, 178, 267, 268 , 269, 270]
append!(final_heights, [[mean(tf[final_index[1:3],:].avg_height), 
                        mean(tf[final_index[4:5],:].avg_height)]])
append!(heights_error, [[std(tf[final_index[1:3],:].avg_height), 
                        std(tf[final_index[4:6],:].avg_height)]])
append!(comparison_time, mean(tf[final_index, :].time))
##
nam = ["A. veronii", "E. coli", "S. cerevisiae (aa)", "S. cerevisiae", "V. cholerae (wt)", "V. cholerae (EPS-)", "K. pneumoniae", "B. cereus", "S. aureus"]
ctg = repeat(["Interferometry", "Control"], inner = 0)
#
groupedbar(nam, reduce(hcat, final_heights)', err=reduce(hcat, heights_error)',
           label = ["Interferometry" "Control"], ylabel="Average Height [μm]", 
           xrotation=40, bottom_margin=3mm, size=(400, 280), grid=false,
           color=[ColorSchemes.okabe_ito[1] :gray], ylim=(0, 520))
#savefig("figs/fig1/interferometry_control.pdf")
