using StatGeochem, Plots

geol = importdataset("geology_combined.csv", importas=:Tuple)
crat = importdataset("craters.csv", importas=:Tuple)

## --- Plot overall histogram as a check
histogram(log10.(crat.DiamKM), # this is log 2 though right?
    framestyle=:box,
    yscale=:log10, 
    xlabel="Log10 crater diameter",
    ylabel="N",
    label="",
)

## -- try calculating an isochron

unit = "AHi"
t = crat.Unit .== unit
step = 0.1
diam_binedges = 0:step:log2(100)
diam_bincenters = (diam_binedges[1:end-1]+diam_binedges[2:end])/2
N = histcounts(log2.(crat.DiamKM[t]), diam_binedges)
density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step

# Plot crater density vs diameter
scatter(2.0.^diam_bincenters, density, 
    framestyle=:box,
    yscale=:log10, 
    xscale=:log10,
    label="",
    ylabel="Crater density [Craters/km^2]",
    xlabel = "Crater diameter [km]",
)

## -- Try plotting isochrons for each unit

for unit in geol.Unit
    t = crat.Unit .== unit
    step = 0.1
    diam_binedges = 0:step:log2(100)
    diam_bincenters = (diam_binedges[1:end-1]+diam_binedges[2:end])/2
    N = histcounts(log2.(crat.DiamKM[t]), diam_binedges)
    density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step

    # Plot crater density vs diameter
    density[density.<=0] .= NaN
    h = scatter(2.0.^diam_bincenters, density, 
        framestyle=:box,
        yscale=:log10, 
        xscale=:log10,
        label="",
        ylabel="Crater density [Craters/km^2]",
        xlabel = "Crater diameter [km]",
    )
    savefig(h, "$unit isochron.pdf")
end
