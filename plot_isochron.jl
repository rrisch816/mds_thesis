using StatGeochem, Plots, Distributions

cd(@__DIR__)
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

## -- Try plotting and fitting isochrons for each unit

# Load craterfreq and mcmc functions
include("CraterModel.jl")

# Tuple to store results
nunits = length(geol.Unit)
ds = (;geol...,
    crater_count = fill(NaN, nunits),
    model_age = fill(NaN, nunits),
    model_age_sigma = fill(NaN, nunits),
    model_erosion = fill(NaN, nunits),
    model_erosion_sigma = fill(NaN, nunits),
    model_erosion_025CI = fill(NaN, nunits),
    model_erosion_975CI = fill(NaN, nunits),
)

step = 0.1
logdiam_binedges = 0:step:log2(1000)
logdiam_bincenters = (logdiam_binedges[1:end-1]+logdiam_binedges[2:end])/2
diam_bincenters = exp2.(logdiam_bincenters)

for i in eachindex(ds.Unit) # loops through each unit
    unit = ds.Unit[i]
    t = crat.Unit .== unit # selects all craters for that specific unit
    ds.crater_count[i] = count(t)
    if count(t) > 1000 # Only consider units with more than 1000 counted craters
        N = histcounts(log2.(crat.DiamKM[t]), logdiam_binedges)
        density = N ./ ds.Area[findfirst(ds.Unit.==unit)] ./ step
        density_sigma = sqrt.(N) ./ ds.Area[findfirst(ds.Unit.==unit)] ./ step 

        # Plot crater density vs diameter
        density[density.<=0] .= NaN # replaces 0 or negative densities with NaN to avoid plotting/log errors
        h = scatter(diam_bincenters, density, 
            yerror = 2*density_sigma,
            framestyle=:box,
            yscale=:log10, 
            xscale=:log10,
            label="Observed densities (N = $(count(t)))",
            ylabel="Crater density [Craters/km^2]",
            xlabel = "Crater diameter [km]",
            title = ds.UnitDesc[i]
        )

        @info "Runing MCMC on unit $unit:"
        ageest = (ds.AgeMin[i]+ds.AgeMax[i])/2
        acceptancedist, lldist, agedist, erosiondist = crater_mcmc(diam_bincenters, density, density_sigma, age=ageest)
        
        ds.model_age[i] = age = mean(agedist)
        ds.model_age_sigma[i] = age_sigma = std(agedist)
        ds.model_erosion[i] = erosion = mean(erosiondist)
        ds.model_erosion_sigma[i] = erosion_sigma = std(erosiondist)
        ds.model_erosion_025CI[i] = nanpctile(erosiondist, 2.5)
        ds.model_erosion_975CI[i] = nanpctile(erosiondist, 97.5)

        diam_cf, density_cf = craterfreq(age, erosion)
        plot!(h, diam_cf, density_cf, 
            label = "age=$(round(age, sigdigits=6)), erosion=$(round(erosion, sigdigits=6))",
            legend = :bottomleft,
            ylims = nanextrema(density_cf),
        )
        savefig(h, "$unit isochron.pdf")
        display(h)

        ha = plot(movmean(acceptancedist,100), ylims=(0,1), 
            ylabel="Acceptance (mean of 100)", 
            title = ds.UnitDesc[i],
            framestyle=:box,
            label="", 
        )
        savefig(ha, "$unit acceptancedist.pdf")

        hl = plot(lldist,
            ylabel="Log likelihood", 
            title = ds.UnitDesc[i],
            framestyle=:box,
            label="", 
        )
        savefig(hl, "$unit lldist.pdf")
    end
end
exportdataset(ds, "results.csv")


## --- Plot resulting erosion rates versus fitted age

    c = resize_colormap(viridis, length(ds.Unit))
    h = scatter(ds.model_age, ds.model_erosion,
        xerror = 2*ds.model_age_sigma,
        yerror = (ds.model_erosion - ds.model_erosion_025CI, ds.model_erosion_975CI - ds.model_erosion),
        framestyle=:box,
        yscale=:log10,
        xlabel = "Model Age [Ga]",
        ylabel = "Erosion [nm/a]",
        label = "",
        color = c,
    )
    for i in eachindex(ds.Unit)
        if !isnan(ds.model_age[i]) && !isnan(ds.model_erosion[i])
            annotate!(h, ds.model_age[i]-2*ds.model_age_sigma[i], ds.model_erosion_975CI[i], text("$(ds.Unit[i])", 8, c[i], :right,))
        end
    end 
    savefig(h, "model_age_vs_erosion.pdf")
    display(h)

## --- Erosion rates versus nominal age

    c = resize_colormap(viridis, length(ds.Unit))
    nominal_age = (ds.AgeMax + ds.AgeMin)/2
    nominal_age_sigma = (ds.AgeMax - ds.AgeMin)/4
    h = scatter(nominal_age, ds.model_erosion,
        xerror = 2*nominal_age_sigma,
        yerror = (ds.model_erosion - ds.model_erosion_025CI, ds.model_erosion_975CI - ds.model_erosion),
        framestyle=:box,
        yscale=:log10,
        xlabel = "Nominal Age [Ga]",
        ylabel = "Erosion [nm/a]",
        label = "",
        color = c,
    )
    for i in eachindex(ds.Unit)
        if !isnan(nominal_age[i]) && !isnan(ds.model_erosion[i])
            annotate!(h, nominal_age[i]-2*nominal_age_sigma[i], ds.model_erosion_975CI[i], text("$(ds.Unit[i])", 8, c[i], :right,))
        end
    end   
    savefig(h, "nominal_age_vs_erosion.pdf")
    display(h)


## --- Modelled and nominal age_sigma

    nominal_age = (ds.AgeMax + ds.AgeMin)/2
    nominal_age_sigma = (ds.AgeMax - ds.AgeMin)/4
    h = scatter(nominal_age, ds.model_age,
        xerror = 2*nominal_age_sigma,
        yerror = 2*ds.model_age_sigma,
        framestyle=:box,
        xlabel = "Nominal Age [Ga]",
        ylabel = "Model Age [Ga]",
        label = "",
        color = c,
    )
    for i in eachindex(ds.Unit)
        if !isnan(nominal_age[i]) && !isnan(ds.model_age[i])
            annotate!(h, nominal_age[i]-2*nominal_age_sigma[i], ds.model_age[i]+2*ds.model_age_sigma[i], text("$(ds.Unit[i])", 8, c[i], :right,))
        end
    end    
    savefig(h, "nominal_vs_model_age.pdf")
    display(h)

## --- Try calculating a particular isochron

# unit = "AHi"
unit = "HNt"
t = crat.Unit .== unit # now extracts only craters in unit AHi (boolean mask)
step = 0.1 # defines bin width in log base 2
logdiam_binedges = 0:step:log2(512) # bin edges from log2(0) to log2(100)
logdiam_bincenters = (logdiam_binedges[1:end-1]+logdiam_binedges[2:end])/2 # computes midpoint of each bin
diam_bincenters = exp2.(logdiam_bincenters)
N = histcounts(log2.(crat.DiamKM[t]), logdiam_binedges) #converts filtered crater diameters to log base 2, counts how many craters in each bin
density = N ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step # normalizes counts by area of unit AHi, finds index where unit equals AHi and grabs area
density_sigma = sqrt.(N) ./ geol.Area[findfirst(geol.Unit.==unit)] ./ step 
# further normalize with bin width step = 0.1
density[density.==0] .= NaN

# Plot crater density vs diameter
h = scatter(diam_bincenters, density,
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

age = 3.6
erosion = 1e1
(CF_crater_D, CF_crater_N) = craterfreq(age, erosion)
plot!(h, CF_crater_D, CF_crater_N, label ="age=$age, erosion=$erosion")
savefig(h, "$unit fitted.pdf")
display(h)

## --- MCMC!

acceptancedist, lldist, agedist, erosiondist = crater_mcmc(diam_bincenters, density, density_sigma)

plot(lldist)
plot(agedist)
plot(erosiondist)

# Get results from stationary distribution
age = mean(agedist)
age_sigma = std(agedist)
erosion = mean(erosiondist)
erosion_sigma = std(erosiondist)

# Plot results!
h = scatter(2.0.^logdiam_bincenters, density, # reverses log2(x) of bincenters
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