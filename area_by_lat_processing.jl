## April 18th, 2025
## Splitting craters by unit and lattitude band

using StatGeochem

abl = importdataset("area_by_lat.csv", importas=:Tuple)
crat = importdataset("craters.csv", importas=:Tuple)

# Need string processing for lat_band column in area_by_lat.csv

# assign craters band numbers? dictionary mapping? -17:-34 = 2 (ex)
for i in length(crat)
    if LAT >=

l1 = 0 .< crat.LAT .< 17


tU = crat.Unit .== "AHi"

t .& tU
 # for loops


 # if count t is less than 1000, break or continue (as opposed to trying to fit the isochron model)

 t = (0 .< crat.LAT .< 17) .& (crat.Unit .== "lHl")


 band_edges = collect(-85.0:17.0:85.0)
 latband_dict = Dict{Tuple{Float64, Float64}, Union{Int, Float64}}()
 
 for i in 1:(length(band_edges)-1)
     pair = (band_edges[i], band_edges[i+1])
     latband_dict[pair] = i
 end

 # Add NaN handling
 latband_dict[(NaN, NaN)] = NaN

unique_units = unique(collect(abl.SIM3292_Global_Geology_Unit))

using DataFrames

groups = Dict{String, DataFrame}()

for i in 1:length(unique_units)
    u = unique_units[i]
    latmin = abl.min_lat[i]
    latmax = abl.max_lat[i]
    pair = (latmin, latmax)

    band = get(latband_dict, pair, NaN)

    if isnan(band)
        continue
    end

    t = (crat.LAT .> latmin) .& (crat.LAT .< latmax) .& (crat.Unit .== u)

    if count(t) < 1000
        continue
    end

    df = DataFrame(Unit = crat.Unit[t], LAT = crat.LAT[t], Band = band, Area = )  # include more fields if needed
    group_name = "$(u)$(band)"
    groups[group_name] = df
end