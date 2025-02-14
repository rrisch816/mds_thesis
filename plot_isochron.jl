using StatGeochem, Plots, Distributions

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

    # diam = 2.0.^diam_bincenters

    # lldist, agedist, erosiondist = crater_mcmc(diam, density, density_sigma)
    # age = mean(agedist)
    # age_sigma = std(agedist)
    # erosion = mean(erosiondist)
    # erosion_sigma = std(erosiondist)

    # diam_cf, density_cf = craterfreq(age, erosion)
    # plot!(h, diam_cf, density_cf, label ="age=$age, erosion=$erosion")

    savefig(h, "$unit isochron.pdf")
end

## -- try calculating a particular isochron

# unit = "AHi"
unit = "HNt"
t = crat.Unit .== unit # now extracts only craters in unit AHi (boolean mask)
step = 0.1 # defines bin width in log base 2
diam_binedges = 0:step:log2(512) # bin edges from log2(0) to log2(100)
diam_bincenters = (diam_binedges[1:end-1]+diam_binedges[2:end])/2 # computes midpoint of each bin
N = histcounts(log2.(crat.DiamKM[t]), diam_binedges) #converts filtered crater diameters to log base 2, counts how many craters in each bin
density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step # normalizes counts by area of unit AHi, finds index where unit equals AHi and grabs area
density_sigma = sqrt.(N) ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step 
# further normalize with bin width step = 0.1
density[density.==0] .= NaN

# Plot crater density vs diameter
h = scatter(2.0.^diam_bincenters, density, # reverses log2(x) of bincenters
    yerror=2*density_sigma,
    ylims=nanextrema(density),
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

age = 3.6
erosion = 1e1
(CF_crater_D, CF_crater_N) = craterfreq(age, erosion)
plot!(h, CF_crater_D, CF_crater_N, label ="age=$age, erosion=$erosion")
savefig(h, "$unit fitted.pdf")
display(h)

## --- MCMC!

diam = 2.0.^diam_bincenters
lldist, agedist, erosiondist = crater_mcmc(diam, density, density_sigma)

plot(lldist)
plot(agedist)
plot(erosiondist)

# Get results from stationary distribution
age = mean(agedist)
age_sigma = std(agedist)
erosion = mean(erosiondist)
erosion_sigma = std(erosiondist)

# Plot results!
h = scatter(2.0.^diam_bincenters, density, # reverses log2(x) of bincenters
    yerror=2*density_sigma,
    ylims=nanextrema(density),
    framestyle=:box,
    yscale=:log10, # log scaling
    xscale=:log10,
    label="",
    ylabel="Crater density [Craters/km^2]",
    xlabel = "Crater diameter [km]",
)
diam_cf, density_cf = craterfreq(age, erosion)
plot!(h, diam_cf, density_cf, label ="age=$age, erosion=$erosion")
savefig(h, "$unit fitted.pdf")
display(h)

## ----

# try messing with bins! Not randomly, needs science

# Try with more units, see what age erosion combos work (always old and much erosion?)
# Try units old medium young, late amazonian too

# Reference curve vs this new hump
# Why does Michael paper curve not have humps

# Mess around with craterstats + IDLE

# 1/31: What do curves actually mean? Sorta fits with massive age, big erosion