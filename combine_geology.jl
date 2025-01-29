1+1
## --- 
using StatGeochem

## loading in dataset as geol in tuple format
geol = importdataset("global_geology.csv", importas=:Tuple)

## from geol dataset, access "Unit" column, finding all distinct Unit values, sorting alphabetically
units =  sort!(unique(geol.Unit))

## finds indices of elements in geol.Unit that match elements in "units", returns array of indices
t = findmatches(units, geol.Unit)
combined = (; ## creates named tuple
    Unit = geol.Unit[t], # extracts geological unit names from geol.Unit using indices from t
    UnitDesc = geol.UnitDesc[t],
    # computes the total area for each unique geological unit
    Area = [sum(geol.SphArea_km[geol.Unit .== u]) for u in units] # creates boolean array checking where geol.Unit equals u
)

exportdataset(combined, "geology_combined.csv")
