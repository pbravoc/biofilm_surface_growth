using DataFrames, CSV 
using Statistics

df = DataFrame(CSV.File("data/sims/bootstrap/all_bootstrap.csv"))
gf = groupby(df, :strain)
percentage_border = 0.1
ql(x) = round(quantile(x, [percentage_border])[1], digits=3) # Quantile low
qh(x) = round(quantile(x, [1-percentage_border][1]), digits=3) # Quantile high
gf_summary = combine(gf, [:α=>ql, :α=>qh,
                          :β=>ql, :β=>qh,
                          :L=>ql, :L=>qh,
                          :h_max=>ql, :h_max=>qh])
print(gf_summary)
