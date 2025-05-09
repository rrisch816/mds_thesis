## April 18th, 2025
## Splitting craters by unit and lattitude band

using StatGeochem

abl = importdataset("area_by_lat.csv", importas=:Tuple)
crat = importdataset("craters.csv", importas=:Tuple)

##
band_edges = collect(-85.0:17.0:85.0)
latband_dict = Dict{Tuple{Float64, Float64}, Union{Int, Float64}}()
 
for i in 1:(length(band_edges)-1)
    latband_dict[(band_edges[i], band_edges[i+1])] = i
end
latband_dict[(NaN, NaN)] = NaN

# Unique combinations of (unit, latmin, latmax)
abl_units = collect(zip(abl.SIM3292_Global_Geology_Unit, abl.min_lat, abl.max_lat))
unique_unit_bands = unique(abl_units)

groups = Dict{String, DataFrame}()

for (u, latmin, latmax) in unique_unit_bands
    pair = (latmin, latmax)
    band = get(latband_dict, pair, NaN)
    if isnan(band)
        continue
    end

    # Find the first matching area in abl for this unit/band
    match = findfirst((abl.SIM3292_Global_Geology_Unit .== u) .& (abl.min_lat .== latmin) .& (abl.max_lat .== latmax))
    if isnothing(match)
        continue
    end
    area = abl.SUM_split_area[match]

    # Filter craters that match this unit and lat band
    t = (crat.Unit .== u) .& (crat.LAT .> latmin) .& (crat.LAT .< latmax)
    if count(t) < 1000
        continue
    end

    df = DataFrame(
        Unit = crat.Unit[t],
        LAT = crat.LAT[t],
        Band = fill(band, count(t)),
        Area = fill(area, count(t)),
        Diameter = crat.DiamKM[t]
    )

    group_name = "$(u)$(band)"
    groups[group_name] = df
end

##--- Starting MCMC process
# Load craterfreq and mcmc functions
include("CraterModel.jl")

group_names = collect(keys(groups))
summary_df = DataFrame(
    group_name = group_names,
    crater_count = [nrow(groups[name]) for name in group_names],
    model_age = fill(NaN, length(group_names)),
    model_age_sigma = fill(NaN, length(group_names)),
    model_erosion = fill(NaN, length(group_names)),
    model_erosion_sigma = fill(NaN, length(group_names)),
    model_erosion_025CI = fill(NaN, length(group_names)),
    model_erosion_975CI = fill(NaN, length(group_names))
)

binwidth = 0.1
logdiam_binedges = 0:binwidth:log2(1000)
logdiam_bincenters = (logdiam_binedges[1:end-1]+logdiam_binedges[2:end])/2
diam_bincenters = exp2.(logdiam_bincenters)

for i in eachindex(summary_df.group_name)
    name = summary_df.group_name[i]
    summary_df.crater_count[i] = nrow(groups[name])
    
    N = histcounts(log2.(crat.DiamKM[t]), logdiam_binedges)
    density = N ./ ds.Area[findfirst(ds.Unit.==unit)] ./ binwidth
    density_sigma = sqrt.(N) ./ ds.Area[findfirst(ds.Unit.==unit)] ./ binwidth 

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