##
using Plots, StatsPlots, Colors
using DataFrames, CSV

df = DataFrame(CSV.File("data/timelapses/database.csv"))
##