## Rebecca Risch
## EARS 13 – Introduction to Computational Methods in Earth Sciences
## Final Project, Thesis Offshoot
## June 2025

using StatGeochem, DataFrames, CSV

three_vars = CSV.read("three_new_vars.csv", DataFrame)
# CSV too large, crashed StatGeochem. Used ChatGPT, which suggested the CSV package

names(three_vars)
for name in names(three_vars)
    println(name)
end

three_vars.mid_lat = ((three_vars.lat_banded_results_csv_min_lat .+ three_vars.lat_banded_results_csv_max_lat) ./ 2)
three_vars.abs_lat = abs.(three_vars.mid_lat)

## CORRELATION: 
using Statistics

# Columns of interest
interest_cols = [:Global_Geology_Intersect1_results_csv_model_erosion,
        :drainage_density,
        :abs_lat,
        :avg_slope_table_MEAN,
        :relief_by_unit_relief]

short_labels = Dict(
    :Global_Geology_Intersect1_results_csv_model_erosion => "Erosion",
    :drainage_density => "Drainage Density",
    :avg_slope_table_MEAN => "Avg Slope",
    :relief_by_unit_relief => "Relief",
    :abs_lat => "Lattitude"
)  
subdf = select(three_vars, interest_cols)
subdf_clean = dropmissing(subdf)

cor_matrix = cor(Matrix(subdf_clean))

short_names = [short_labels[c] for c in interest_cols]
cor_df = DataFrame(cor_matrix, short_names)
cor_df.variable = short_names
cor_df = select(cor_df, ["variable"; short_names]) # reordering

for col in short_names
    cor_df[!, col] = round.(cor_df[!, col], digits=2)
end

CSV.write("cor_matrix.csv", cor_df)


## DISTRIBUTIONS of variables of interest:      
subdf = dropmissing(select(three_vars, interest_cols))
summary_df = describe(subdf, :min, :q25, :median, :mean, :q75, :max)
summary_df.variable = [short_labels[Symbol(v)] for v in summary_df.variable]

for col in names(summary_df)
    if eltype(summary_df[!, col]) <: Union{Nothing, Missing, Real} # got this line from ChatGPT, debugging when rounding wasn't working
        summary_df[!, col] = map(x -> x === nothing ? missing : round(x, digits=2), summary_df[!, col])
    end
end

CSV.write("summary_table.csv", summary_df)

## MODEL: 
using CategoricalArrays

using GLM
model = glm(@formula(Global_Geology_Intersect1_results_csv_model_erosion ~ drainage_density + avg_slope_table_MEAN + relief_by_unit_relief), three_vars, Normal(), IdentityLink())
# use to have "+ lat_banded_results_csv_lat_band" in the formula, but that needs banded erosion as its DV

model_table = coeftable(model)
model_df = DataFrame(model_table)

numeric_cols = names(model_df, Not("Name"))

for col in numeric_cols
    if eltype(model_df[!, col]) <: Union{Missing, Real}
        model_df[!, col] = round.(model_df[!, col]; digits=3)
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

lat_df = dropmissing(select(three_vars, [
    :lat_banded_results_csv_model_erosion,
    :abs_lat,
    :Global_Geology_Intersect1_SIM3292_Global_Geology_Unit
]))

# Fit the GLM
lat_model = glm(@formula(lat_banded_results_csv_model_erosion ~ 
                     abs_lat + 
                     Global_Geology_Intersect1_SIM3292_Global_Geology_Unit),
            lat_df, Normal(), IdentityLink())

lat_model_table = coeftable(lat_model)
lat_model_df = DataFrame(lat_model_table)

lat_numerics = names(lat_model_df, Not("Name"))
for col in lat_numerics
    if eltype(lat_model_df[!, col]) <: Union{Missing, Real}
        lat_model_df[!, col] = round.(lat_model_df[!, col]; digits=3)
    end
end

CSV.write("lat_glm_coefficients.csv", lat_model_df)


## Distribution visualizations:
using Plots, StatsPlots, StatsBase

for (col, label) in short_labels
    data = subdf[!, col]

    min, max = minimum(data), maximum(data)
    binedges = range(min, max; length=25)

    h = fit(Histogram, data, binedges)
    x = midpoints(h.edges[1])
    widths = step(binedges)
    y = h.weights ./ (sum(h.weights) * widths)

    p1 = plot(x, y, seriestype=:bar, xlabel=label, ylabel="Density",
              title="PDF of "*label, label="", legend=false)
    savefig(p1, string(label)*"_pdf.png")
    display(p1)

    # LOG SPACE
    logdata = log10.(data)
    min_log, max_log = minimum(logdata), maximum(logdata)
    binedges_log = range(min_log, max_log; length=25)

    h_log = fit(Histogram, logdata, binedges_log)
    x_log = midpoints(h_log.edges[1])
    widths_log = step(binedges_log)
    y_log = h_log.weights ./ (sum(h_log.weights) * widths_log)

    p2 = plot(x_log, y_log, seriestype=:bar, xlabel="log("*label*")", ylabel="Density",
              title="Log10-Scaled PDF of "*label, label="", legend=false)
    savefig(p2, string(label)*"_log_pdf.png")
    display(p2)
end

## SCATTERPLOTS:
using StatsPlots

lat_sub = select(three_vars, :lat_banded_results_csv_model_erosion, :abs_lat)
lat_sub = dropmissing(lat_sub)

lat_plot = @df lat_sub StatsPlots.scatter(
    :abs_lat,
    :lat_banded_results_csv_model_erosion,
    xlabel = "abs(Median Latitude of Latitude Band (ranges 17º each))",
    ylabel = "Erosion (nm/a, log scale)",
    label = "Observed",
    title = "Erosion by Latitude Band (log scale)",
    yscale = :log10,
    ylims = (10e-2, 10e3),
    legend = true
)
savefig(lat_plot, "lat_plot.png")

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

dd_sub.log_erosion = log10.(dd_sub.Global_Geology_Intersect1_results_csv_model_erosion)
dd_log_scatter = @df dd_sub Plots.scatter(
    :drainage_density,
    :log_erosion,
    xlabel = "Drainage Density (km/km^2)",
    ylabel = "Log Base 10 Erosion Rate log(nm/a)",
    label = "Observed",
    legend = true,
    title = "Log Erosion vs Drainage Density",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3,
)
savefig(dd_log_scatter, "dd_log_scatter.png")

relief_sub = select(three_vars, :Global_Geology_Intersect1_results_csv_model_erosion, :relief_by_unit_relief)
relief_sub = dropmissing(relief_sub)

rl_scatter = @df relief_sub Plots.scatter(
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
savefig(rl_scatter, "relief_scatter.png")

relief_sub.log_erosion = log10.(relief_sub.Global_Geology_Intersect1_results_csv_model_erosion)
rl_log_scatter = @df relief_sub Plots.scatter(
    :relief_by_unit_relief,
    :log_erosion,
    xlabel = "Relief (m)",
    ylabel = "Log Base 10 Erosion Rate log(nm/a)",
    label = "Observed",
    legend = true,
    title = "Log Erosion vs Relief",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3,
)
savefig(rl_log_scatter, "relief_log_scatter.png")

slope_sub = select(three_vars, :Global_Geology_Intersect1_results_csv_model_erosion, :avg_slope_table_MEAN)
slope_sub = dropmissing(slope_sub)

sl_scatter = @df slope_sub Plots.scatter(
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
savefig(sl_scatter, "slope_scatter.png")

slope_sub.log_erosion = log10.(slope_sub.Global_Geology_Intersect1_results_csv_model_erosion)
sl_log_scatter = @df slope_sub Plots.scatter(
    :avg_slope_table_MEAN,
    :log_erosion,
    xlabel = "Slope (º)",
    ylabel = "Log Base 10 Erosion Rate log(nm/a)",
    label = "Observed",
    legend = true,
    title = "Log Erosion vs Slope",
    jitter = 0.3,
    alpha = 0.4,
    markersize = 3,
)
savefig(sl_log_scatter, "slope_log_scatter.png")


## PREDICTING RESULTS
## Lat band model results viz:

lat_sub.predicted_erosion = predict(lat_model)

lat_results_overlay = @df lat_sub StatsPlots.scatter!(lat_plot,
    :abs_lat,
    :predicted_erosion,
    label = "Predicted",
    color = :pink,
    yscale = :log10,
    ylims = (10e-2, 10e3),
    legend = true
)

savefig(lat_results_overlay, "lat_plot_pred.png")

## Drainage Density Model Results Viz:

xrange = range(
    minimum(three_vars.drainage_density),
    stop = maximum(three_vars.drainage_density),
    length = 100
)

mean_sl = mean(skipmissing(three_vars.avg_slope_table_MEAN))
mean_rl = mean(skipmissing(three_vars.relief_by_unit_relief))

dd_pred_df = DataFrame(
    drainage_density = collect(xrange),
    avg_slope_table_MEAN = fill(mean_sl, 100),
    relief_by_unit_relief = fill(mean_rl, 100)
)

dd_pred_df.predicted_erosion = predict(model, dd_pred_df)

dd_modeled_results = @df dd_pred_df StatsPlots.plot!(dd_scatter,
    :drainage_density,
    :predicted_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(dd_modeled_results, "dd_modeled_results.png")

dd_pred_df.log_pred_erosion = log10.(dd_pred_df.predicted_erosion)

dd_log_modeled_results = @df dd_pred_df StatsPlots.plot!(dd_log_scatter,
    :drainage_density,
    :log_pred_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(dd_log_modeled_results, "dd_log_modeled_results.png")

## Relief Model Results Viz:

rlxrange = range(
    minimum(three_vars.relief_by_unit_relief),
    stop = maximum(three_vars.relief_by_unit_relief),
    length = 100
)

mean_dd = mean(skipmissing(three_vars.drainage_density))

rl_pred_df = DataFrame(
    drainage_density = fill(mean_dd, 100),
    avg_slope_table_MEAN = fill(mean_sl, 100),
    relief_by_unit_relief = collect(rlxrange)
)

rl_pred_df.predicted_erosion = predict(model, rl_pred_df)

rl_modeled_results = @df rl_pred_df StatsPlots.plot!(rl_scatter,
    :relief_by_unit_relief,
    :predicted_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(rl_modeled_results, "rl_modeled_results.png")

rl_pred_df.log_pred_erosion = log10.(rl_pred_df.predicted_erosion)

rl_log_modeled_results = @df rl_pred_df StatsPlots.plot!(rl_log_scatter,
    :relief_by_unit_relief,
    :log_pred_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(rl_log_modeled_results, "rl_log_modeled_results.png")

# test_rl_model = glm(@formula(Global_Geology_Intersect1_results_csv_model_erosion ~ relief_by_unit_relief), three_vars, Normal(), IdentityLink())
# rl_pred_df2 = DataFrame(
#     relief_by_unit_relief = collect(rlxrange)
# )
# rl_pred_df2.predicted_erosion2 = predict(test_rl_model, rl_pred_df2)
# rl_modeled_results = @df rl_pred_df2 StatsPlots.scatter!(rl_scatter,
#     :relief_by_unit_relief,
#     :predicted_erosion2,
#     label = "Predicted",
#     color = :red,
#     legend = true
# )

## Slope Model Results Viz:

slxrange = range(
    minimum(three_vars.avg_slope_table_MEAN),
    stop = maximum(three_vars.avg_slope_table_MEAN),
    length = 100
)

sl_pred_df = DataFrame(
    drainage_density = fill(mean_dd, 100),
    avg_slope_table_MEAN = collect(slxrange),
    relief_by_unit_relief = fill(mean_sl, 100)
)

sl_pred_df.predicted_erosion = predict(model, sl_pred_df)

sl_modeled_results = @df sl_pred_df StatsPlots.plot!(sl_scatter,
    :avg_slope_table_MEAN,
    :predicted_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(sl_modeled_results, "sl_modeled_results.png")

sl_pred_df.log_pred_erosion = log10.(sl_pred_df.predicted_erosion)

sl_log_modeled_results = @df sl_pred_df StatsPlots.plot!(sl_log_scatter,
    :avg_slope_table_MEAN,
    :log_pred_erosion,
    label = "Predicted",
    color = :red,
    legend = true
)

savefig(sl_log_modeled_results, "sl_log_modeled_results.png")

## MODEL EVALUATION:
model_eval_df = dropmissing(select(three_vars, [
    :Global_Geology_Intersect1_results_csv_model_erosion,
    :drainage_density,
    :avg_slope_table_MEAN,
    :relief_by_unit_relief
]))

y_actual = model_eval_df.Global_Geology_Intersect1_results_csv_model_erosion
y_pred = predict(model, model_eval_df)
resid = y_actual .- y_pred

ss_total = sum((y_actual .- mean(y_actual)).^2)
ss_res = sum(resid .^ 2)
r_squared = 1 - ss_res / ss_total

rmse = sqrt(mean(resid .^ 2))

println("R² = ", round(r_squared, digits=4))
println("RMSE = ", round(rmse, digits=2))

## LOG MODEL AND EVAL
three_vars.log_erosion = log10.(three_vars.Global_Geology_Intersect1_results_csv_model_erosion)

log_model = glm(@formula(log_erosion ~ drainage_density + avg_slope_table_MEAN + relief_by_unit_relief),
                three_vars, Normal(), IdentityLink())

log_model_eval_df = dropmissing(select(three_vars, [
    :log_erosion,
    :drainage_density,
    :avg_slope_table_MEAN,
    :relief_by_unit_relief
]))

y_actual = log_model_eval_df.log_erosion
y_pred = predict(log_model, log_model_eval_df)
resid = y_actual .- y_pred

ss_total = sum((y_actual .- mean(y_actual)).^2)
ss_res = sum(resid .^ 2)
r_squared = 1 - ss_res / ss_total

rmse = sqrt(mean(resid .^ 2))

println("R² = ", round(r_squared, digits=4))
println("RMSE = ", round(rmse, digits=2))

## NEW LAT MODELS:

new_lat_df = dropmissing(select(three_vars, [
    :lat_banded_results_csv_model_erosion,
    :abs_lat,
    :Global_Geology_Intersect1_SIM3292_Global_Geology_Unit,
    :drainage_density
]))

lat_model2 = glm(@formula(lat_banded_results_csv_model_erosion ~ 
                     abs_lat + drainage_density),
                     new_lat_df, Normal(), IdentityLink())

lat_model2_table = coeftable(lat_model2)
lat_model2_df = DataFrame(lat_model2_table)

for col in names(lat_model2_df)
    if eltype(lat_model2_df[!, col]) <: Real
        lat_model2_df[!, col] = string.(round.(lat_model2_df[!, col], sigdigits=4))
    end
end
CSV.write("lat_glm_coefficients2.csv", lat_model2_df)

lat_model_eval = dropmissing(select(three_vars, [
    :lat_banded_results_csv_model_erosion,
    :abs_lat,
    :drainage_density
]))

y_actual = lat_model_eval.lat_banded_results_csv_model_erosion
y_pred = predict(lat_model2, lat_model_eval)
resid = y_actual .- y_pred

ss_total = sum((y_actual .- mean(y_actual)).^2)
ss_res = sum(resid .^ 2)
r_squared = 1 - ss_res / ss_total
rmse = sqrt(mean(resid .^ 2))

println("R² = ", round(r_squared, digits=4))
println("RMSE = ", round(rmse, digits=2))


## DD vs Latitude
grouped = combine(groupby(new_lat_df, :abs_lat),
    nrow => :Count,
    :drainage_density => mean => :mean_dd,
    :drainage_density => std => :sd_dd)
show(grouped, allrows=true)

lat_model3 = glm(@formula(drainage_density ~ abs_lat), new_lat_df, Normal(), IdentityLink())

lat_model3_table = coeftable(lat_model3)
lat_model3_df = DataFrame(lat_model3_table)

for col in names(lat_model3_df)
    if eltype(lat_model3_df[!, col]) <: Real
        lat_model3_df[!, col] = string.(round.(lat_model3_df[!, col], sigdigits=4))
    end
end
CSV.write("lat_glm_coefficients3.csv", lat_model3_df)
    

lat_dd_model_eval = dropmissing(select(three_vars, [
    :abs_lat,
    :drainage_density
]))

y_actual = lat_dd_model_eval.drainage_density
y_pred = predict(lat_model3, lat_dd_model_eval)
resid = y_actual .- y_pred

ss_total = sum((y_actual .- mean(y_actual)).^2)
ss_res = sum(resid .^ 2)
r_squared = 1 - ss_res / ss_total
rmse = sqrt(mean(resid .^ 2))

println("R² = ", round(r_squared, digits=4))
println("RMSE = ", round(rmse, digits=2))

## WILNER DATA SET STUFF
using StatGeochem
glacial = importdataset("glacial_erosion_Mars.tsv", '\t', importas=:Tuple)
non = importdataset("nonglacial_erosion_Mars.tsv", '\t', importas=:Tuple)

# Log data that brenhin gave- lots of unit changes
# log time scale vs log rate
# presence or absence of saddlerian trends
# Negative trends typically in earth processes, though not glacial erosion
# Convert 
# Rate is in mm/yr for rate, years for time interval. Convert from billions of years to years, nm to mm

using DataFrames

glacial_df = DataFrame(glacial)
glacial_df.erosion_converted = 1_000_000 .* glacial_df.Erosion_rate_mm_yr
glacial_df.log_erosion = log10.(glacial_df.erosion_converted)
glacial_df.log_time = log10.(glacial_df.Time_interval_yr)

non_df = DataFrame(non)
non_df.erosion_converted = 1_000_000 .* non_df.Erosion_rate_mm_yr
non_df.log_erosion = log10.(non_df.erosion_converted)
non_df.log_time = log10.(non_df.Time_interval_yr)

glacial_plot = @df glacial_df Plots.scatter(
    :log_time,
    :log_erosion,
    xlabel = "Log Base 10 Time Scale (yr)",
    ylabel = "Log Base 10 Erosion rate (nm/yr)",
    legend = true,
    title = "Time Scale vs Erosion Rate (Log Base 10)",
    jitter = 0.3,
    alpha = 0.8,
    label = "Glacial",
    markersize = 3
)

non_plot = @df non_df Plots.scatter!(glacial_plot,
    :log_time,
    :log_erosion,
    legend = true,
    jitter = 0.3,
    alpha = 0.8,
    markersize = 3,
    col = "red",
    label = "Non Glacial"
)

savefig(non_plot, "sadler_wilner.png")

three_vars.log_time = log10.(three_vars.Global_Geology_Intersect1_results_csv_model_age .* 1_000_000_000)
three_vars.log_erosion = log10.(three_vars.Global_Geology_Intersect1_results_csv_model_erosion)
three_vars.log_banded_erosion = log10.(three_vars.lat_banded_results_csv_model_erosion)

my_plot = scatter!(non_plot,
    three_vars.log_time,
    three_vars.log_erosion,
    legend = true,
    jitter = 0.3,
    alpha = 0.8,
    markersize = 3,
    color = "seagreen4",
    label = "Risch Data"
)

savefig(my_plot, "my_plot.png")

my_iso_plot = scatter(three_vars.log_time,
    three_vars.log_erosion,
    legend = true,
    jitter = 0.3,
    alpha = 0.8,
    markersize = 4,
    color = "seagreen4",
    label = "Risch Data",
    xlabel = "Log Base 10 Time Scale (yr)",
    ylabel = "Log Base 10 Erosion rate (nm/yr)",
    title = "Time Scale vs Erosion Rate (Log Base 10)",
)

savefig(my_iso_plot, "my_iso_plot.png")

cleaned = dropmissing(three_vars, :abs_lat) # ChatGPT debugging suggestion

my_iso_lat_plot = scatter(cleaned.log_time,
    cleaned.log_banded_erosion,
    legend = false,
    jitter = 0.3,
    alpha = 0.8,
    markersize = 4,
    zcolor = cleaned.abs_lat,
    colormap = :viridis,
    colorbar = true,
    xlabel = "Log Base 10 Time Scale (yr)",
    ylabel = "Log Base 10 Latitude Banded Erosion rate (nm/yr)",
    title = "Time Scale vs Erosion Rate (Log Base 10)",
    guidefont = font(10)
)

savefig(my_iso_lat_plot, "my_iso_lat_plot.png")

# INTERACTION MODEL
interaction_model = glm(@formula(Global_Geology_Intersect1_results_csv_model_erosion ~ drainage_density + avg_slope_table_MEAN + relief_by_unit_relief + drainage_density * avg_slope_table_MEAN + drainage_density * relief_by_unit_relief + avg_slope_table_MEAN * relief_by_unit_relief), three_vars, Normal(), IdentityLink())

int_table = coeftable(interaction_model)
int_df = DataFrame(int_table)

numeric_cols = names(int_df, Not("Name"))

for col in numeric_cols
    if eltype(int_df[!, col]) <: Union{Missing, Real}
        int_df[!, col] = round.(int_df[!, col]; digits=3)
    end
end

CSV.write("int_glm_coefficients.csv", int_df)

int_eval_df = dropmissing(select(three_vars, [
    :Global_Geology_Intersect1_results_csv_model_erosion,
    :drainage_density,
    :avg_slope_table_MEAN,
    :relief_by_unit_relief
]))

y_actual = int_eval_df.Global_Geology_Intersect1_results_csv_model_erosion
y_pred = predict(interaction_model, int_eval_df)
resid = y_actual .- y_pred

ss_total = sum((y_actual .- mean(y_actual)).^2)
ss_res = sum(resid .^ 2)
r_squared = 1 - ss_res / ss_total

rmse = sqrt(mean(resid .^ 2))

println("R² = ", round(r_squared, digits=4))
println("RMSE = ", round(rmse, digits=2))

# martian evidence of sadler effect- discussion, refer wilner et al

# after accounting for sadler effect, glacial erosion rates still higher?


# add in spread summary

# in report, glacial non glacial data, how was in compiled (see wilner et al paper: "Materials and Methods")


## Plot my data over the sadler graph, measured age vs measured erosion rate both in log base 10

## Plot the points erosion by latitude also, color by lat band


## Absolute value of latitude, median (midpoint) of lat band. On continuous scale, not categorical
# treat 90 as -90, abs value

## Two sadler plots (whole erosion and sliced by lat erosion)
## Sliced by lat vs 3 geomorph vars

## normalize and standardize?
## for comparability, standardize vars (- mean / std dev, all on same scale)
## then bigger coeffs means matters more

## dont need to refer to all figures in text, some can be supplementary
