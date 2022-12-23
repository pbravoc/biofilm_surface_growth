#=
This is the figure that shows all the parameters, h_max 
for all the strains in the same plot. Maybe there's a way to 
join this figure w number 6 (as insets?)
=#
##
using DataFrames, CSV
using Statistics
## Load bootstrap for statistics 
Df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
Df = filter(x->x.α .> 0 && x.β .>0 && x.L .>0, Df)                  
gf = groupby(Df, :strain)                           # Group by strain
##
percentage_border = 0.025                           # 2.5% on top + 2.5% on bottom = 5%
ql(x) = round(quantile(x, [percentage_border])[1], digits=3)    # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3)  # Quantile high
df = combine(gf, [:α=>ql, :α=>qh, :α=>mean,
                  :β=>ql, :β=>qh, :β=>mean,
                  :L=>ql, :L=>qh, :L=>mean,
                  :h_max=>ql, :h_max=>qh])
print(df)
df.α_mean = round.(df.α_mean, digits=3)
df.β_mean = round.(df.β_mean, digits=3)
df.L_mean = round.(df.L_mean, digits=3)

## Load best-fit parameters 
Df = DataFrame(CSV.File("data/timelapses/fit_params_interface.csv"))
Df = filter(x->x.fit .== "48h", Df)
df[:, :α] = round.(Df.x1, digits=3)
df[:, :β] = round.(Df.x2, digits=3)
df[:, :L] = round.(Df.x3, digits=3)

## Add species 
df[:, :name] = ["A. veronii", "E. coli", "S. cerevisiae (aa)",
                "S. cerevisiae (wt)", "V. cholerae (wt)", "V. cholerae (EPS-)",
                "K. pneumoniae", "B. cereus", "S. aureus"]

# Reorder for simplicity
df = df[!, [16, 13, 14, 15, 2, 3, 5, 6, 8, 9, 4, 7, 10, 1]]
CSV.write("data/table1_data.csv", df)