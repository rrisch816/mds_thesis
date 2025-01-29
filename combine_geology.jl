1+1
## --- 
using StatGeochem

geol = importdataset("global_geology.csv", importas=:Tuple)

units =  sort!(unique(geol.Unit))
t = findmatches(units, geol.Unit)
combined = (;
    Unit = geol.Unit[t],
    UnitDesc = geol.UnitDesc[t],
    Area = [sum(geol.SphArea_km[geol.Unit .== u]) for u in units] # geol.SUM_SphArea_km[t],
)

exportdataset(combined, "geology_combined.csv")
