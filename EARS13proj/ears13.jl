## Rebecca Risch
## EARS 13 – Introduction to Computational Methods in Earth Sciences
## Final Project, Thesis Offshoot
## June 2025

using StatGeochem, DataFrames, CSV

three_vars = CSV.read("three_new_vars.csv", DataFrame)

names(three_vars)
for name in names(three_vars)
    println(name)
end


## CORRELATION: 
using Statistics

# Columns of interest
interest_cols = [:Global_Geology_Intersect1_results_csv_model_erosion,
        :drainage_density,
        :lat_banded_results_csv_lat_band,
        :avg_slope_table_MEAN,
        :relief_by_unit_relief]

short_labels = Dict(
    :Global_Geology_Intersect1_results_csv_model_erosion => "Erosion",
    :drainage_density => "Drainage Density",
    :avg_slope_table_MEAN => "Avg Slope",
    :relief_by_unit_relief => "Relief",
    :lat_banded_results_csv_lat_band => "Lat Band"
)  
subdf = select(three_vars, interest_cols)
subdf_clean = dropmissing(subdf)

cor_matrix = cor(Matrix(subdf_clean))

short_names = [short_labels[c] for c in interest_cols]
cor_df = DataFrame(cor_matrix, short_names)
cor_df.variable = short_names
cor_df = select(cor_df, ["variable"; short_names])

for col in short_names
    cor_df[!, col] = round.(cor_df[!, col], digits=2)
end

CSV.write("cor_matrix.csv", cor_df)


## DISTRIBUTIONS of variables of interest:      
subdf = dropmissing(select(three_vars, interest_cols))
summary_df = describe(subdf, :min, :q25, :median, :mean, :q75, :max)
summary_df.variable = [short_labels[Symbol(v)] for v in summary_df.variable]

for col in names(summary_df)
    if eltype(summary_df[!, col]) <: Union{Nothing, Missing, Real}
        summary_df[!, col] = map(x -> x === nothing ? missing : round(x, digits=2), summary_df[!, col])
    end
end

CSV.write("summary_table.csv", summary_df)

## MODEL: 
using CategoricalArrays
# Convert lat_band and unit to categorical
three_vars.lat_banded_results_csv_lat_band = categorical(three_vars.lat_banded_results_csv_lat_band)

using GLM
model = glm(@formula(Global_Geology_Intersect1_results_csv_model_erosion ~ drainage_density + avg_slope_table_MEAN + relief_by_unit_relief), three_vars, Normal(), IdentityLink())
# use to have "+ lat_banded_results_csv_lat_band" in the formula, but that needs banded erosion as its DV

model_table = coeftable(model)
model_df = DataFrame(model_table)

for col in names(model_df)
    if eltype(model_df[!, col]) <: Real
        model_df[!, col] = string.(round.(model_df[!, col], sigdigits=4))
    end
end
CSV.write("whole_glm_coefficients.csv", model_df)


## Average slope is by:
sorted = sort(three_vars, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit,  :lat_banded_results_csv_lat_band, :avg_slope_table_MEAN])

subset = select(sorted, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit, :avg_slope_table_MEAN, :lat_banded_results_csv_lat_band])
show(first(subset, 300), allrows=true, allcols=true)
# UNIT

## Average relief is by:
sorted = sort(three_vars, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit,  :lat_banded_results_csv_lat_band, :relief_by_unit_relief])

subset = select(sorted, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit, :relief_by_unit_relief, :lat_banded_results_csv_lat_band])
show(first(subset, 300), allrows=true, allcols=true)
# UNIT

# Drainage density is by:
sorted = sort(three_vars, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit,  :lat_banded_results_csv_lat_band, :drainage_density])

subset = select(sorted, [:Global_Geology_Intersect1_SIM3292_Global_Geology_Unit, :drainage_density, :lat_banded_results_csv_lat_band])
show(first(subset, 300), allrows=true, allcols=true)
# UNIT
# So for these, we use whole erosion results (not banded unit combo results)

# LAT BAND needs its own model then
three_vars.Global_Geology_Intersect1_SIM3292_Global_Geology_Unit = categorical(three_vars.Global_Geology_Intersect1_SIM3292_Global_Geology_Unit)

# Drop missing values in the relevant columns
df = dropmissing(select(three_vars, [
    :lat_banded_results_csv_model_erosion,
    :lat_banded_results_csv_lat_band,
    :Global_Geology_Intersect1_SIM3292_Global_Geology_Unit
]))

# Fit the GLM
lat_model = glm(@formula(lat_banded_results_csv_model_erosion ~ 
                     lat_banded_results_csv_lat_band + 
                     Global_Geology_Intersect1_SIM3292_Global_Geology_Unit),
            df, Normal(), IdentityLink())

lat_model_table = coeftable(lat_model)
lat_model_df = DataFrame(lat_model_table)

for col in names(lat_model_df)
    if eltype(lat_model_df[!, col]) <: Real
        lat_model_df[!, col] = string.(round.(lat_model_df[!, col], sigdigits=4))
    end
end
CSV.write("lat_glm_coefficients.csv", lat_model_df)



## Distribution visualizations:
using Plots, StatsPlots, StatsBase

interest_cols = keys(short_labels)  # ensures only those columns are used

for col in interest_cols
    data = subdf[!, col]

    minval, maxval = minimum(data), maximum(data)
    binedges = range(minval, maxval; length=25)

    h = fit(Histogram, data, binedges)
    x = midpoints(h.edges[1])
    y = h.weights ./ (sum(h.weights) * step(binedges))  # normalize to PDF

    # Use the short label for plotting
    label = short_labels[col]
    p = plot(x, y, seriestype=:bar, xlabel=label, ylabel="Density",
             title="PDF of "*label, label="", legend=false)

    savefig(p, string(label)*"_pdf.png")
    display(p)
end

## SCATTERPLOTS:
using StatsPlots

latband_sub = select(three_vars, :lat_banded_results_csv_model_erosion, :lat_banded_results_csv_lat_band)
latband_sub = dropmissing(latband_sub)

latband_boxplot = @df latband_sub StatsPlots.boxplot(
    :lat_banded_results_csv_lat_band,
    :lat_banded_results_csv_model_erosion,
    xlabel = "Latitude Band (ranges 17º each)",
    ylabel = "Erosion (nm/a, log scale)",
    label = "Observed",
    title = "Erosion by Latitude Band (log scale)",
    yscale = :log10,
    ylims = (10e-2, 10e3),
    legend = true
)
savefig(latband_boxplot, "latband_boxplot.png")

using Plots

dd_sub = select(three_vars, :Global_Geology_Intersect1_results_csv_model_erosion, :drainage_density)
dd_sub = dropmissing(dd_sub)

dd_scatter = @df dd_sub Plots.scatter(
    :drainage_density,
    :Global_Geology_Intersect1_results_csv_model_erosion,
    xlabel = "Drainage Density (km/km^2)",
    ylabel = "Erosion (nm/a)",
    label = "Observed",
    legend = true,
    title = "Erosion vs Drainage Density",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3
)
savefig(dd_scatter, "dd_scatter.png")

relief_sub = select(three_vars, :Global_Geology_Intersect1_results_csv_model_erosion, :relief_by_unit_relief)
relief_sub = dropmissing(relief_sub)

@df relief_sub Plots.scatter(
    :relief_by_unit_relief,
    :Global_Geology_Intersect1_results_csv_model_erosion,
    xlabel = "Relief (m)",
    ylabel = "Erosion (nm/a)",
    legend = false,
    title = "Erosion vs Relief",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3
)
savefig("relief_scatter.png")

slope_sub = select(three_vars, :Global_Geology_Intersect1_results_csv_model_erosion, :avg_slope_table_MEAN)
slope_sub = dropmissing(slope_sub)

@df slope_sub Plots.scatter(
    :avg_slope_table_MEAN,
    :Global_Geology_Intersect1_results_csv_model_erosion,
    xlabel = "Slope (º)",
    ylabel = "Erosion (nm/a)",
    legend = false,
    title = "Erosion vs Slope",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3
)
savefig("slope_scatter.png")

## Lat band model results viz:

df.predicted_erosion = predict(lat_model)

latband_results_overlay = @df df StatsPlots.boxplot!(latband_boxplot,
    :lat_banded_results_csv_lat_band,
    :predicted_erosion,
    label = "Predicted",
    fillalpha = 0.3,
    linecolor = :red,
    mediancolor = :red,
    outliercolor = :red,
    color = :red,
    yscale = :log10,
    ylims = (10e-2, 10e3),
    legend = true
)
ENV["GKSwstype"] = "100"
savefig(latband_results_overlay, "latband_box_pred.pdf")

## Drainage Density Model Results Viz:

xrange = range(
    minimum(three_vars.drainage_density),
    stop = maximum(three_vars.drainage_density),
    length = 100
)

mean_slope = mean(skipmissing(three_vars.avg_slope_table_MEAN))
mean_relief = mean(skipmissing(three_vars.relief_by_unit_relief))

pred_df = DataFrame(
    drainage_density = collect(xrange),
    avg_slope_table_MEAN = fill(mean_slope, 100),
    relief_by_unit_relief = fill(mean_relief, 100)
)

pred_df.predicted_erosion = predict(model, pred_df)

dd_modeled_results = @df pred_df StatsPlots.scatter!(dd_scatter,
    :drainage_density,
    :predicted_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(dd_modeled_results, "dd_modeled_results.png")

# NEXT: Result interpretation, graphs/viz!!
# Plot data points erosion rate vs lat, erosion rate vs dd, etc 2D scatterplots
# and plot GLMS overlying those?

# Log data that brenhin gave- lots of unit changes
# log time sccale vs log rate
# presence or absence of saddlerian trends
# Negative trends typically in earth processes, though not glacial erosion
# Convert 
# Rate is in mm/yr for rate, years for time interval. Convert from billions of years to years, nm to mm