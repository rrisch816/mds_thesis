using StatGeochem, Plots

geol = importdataset("geology_combined.csv", importas=:Tuple)
crat = importdataset("craters.csv", importas=:Tuple)

## --- Plot overall histogram as a check
histogram(log10.(crat.DiamKM), # creates a histogram of crater diameters on log scale (base 10)
    framestyle=:box, # encloses axes
    yscale=:log10, # sets y axis to log 10 scale
    xlabel="Log10 crater diameter",
    ylabel="N",
    label="", # removes legend label
)

## -- Try plotting isochrons for each unit

for unit in geol.Unit # loops through each unit
    t = crat.Unit .== unit # selects all craters for that specific unit
    step = 0.1
    diam_binedges = 0:step:log2(1000)
    diam_bincenters = (diam_binedges[1:end-1]+diam_binedges[2:end])/2
    N = histcounts(log2.(crat.DiamKM[t]), diam_binedges)
    density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step

    # Plot crater density vs diameter
    density[density.<=0] .= NaN # replaces 0 or negative densities with NaN to avoid plotting/log errors
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

## -- try calculating a particular isochron

# unit = "AHi"
unit = "mNh"
t = crat.Unit .== unit # now extracts only craters in unit AHi (boolean mask)
step = 0.1 # defines bin width in log base 2
diam_binedges = 0:step:log2(512) # bin edges from log2(0) to log2(100)
diam_bincenters = (diam_binedges[1:end-1]+diam_binedges[2:end])/2 # computes midpoint of each bin
N = histcounts(log2.(crat.DiamKM[t]), diam_binedges) #converts filtered crater diameters to log base 2, counts how many craters in each bin
density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step # normalizes counts by area of unit AHi, finds index where unit equals AHi and grabs area
# further normalize with bin width step = 0.1
density[density.==0] .= NaN

# Plot crater density vs diameter
h = scatter(2.0.^diam_bincenters, density, # reverses log2(x) of bincenters
    framestyle=:box,
    yscale=:log10, # log scaling
    xscale=:log10,
    label="",
    ylabel="Crater density [Craters/km^2]",
    xlabel = "Crater diameter [km]",
)

## --- Plot this against a calulated isochron

# To get craterfreq function
include("CraterModel.jl")

(CF_crater_D, CF_crater_N) = craterfreq(4.5, 0.02, 512, true, false, 1e4)
plot!(h, CF_crater_D, CF_crater_N)