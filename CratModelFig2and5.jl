## --- Load required packages

    using Statistics, StatsBase, Plots
    using ProgressMeter: @showprogress

## --- Define craterfreq function

"""
'''julia
function craterfreq(age::Number=1, diameter_min::Number=0.0039,
    \tdiameter_max::Number=1024, EROSION::Bool=false,
    \tSECONDARY::Bool=false, Beta::Number=1e-5; HARTMANN_PROD::Bool=true)
```

Based on range of craters chosen and surface age (Ga), returns list of crater
frequencies for 1 Ga surface age (no. craters/km2). If EROSION is true uses
Smith et al 2008 model for erosion and B is erosion rate (nm/a).

If HARTMANN_PROD if false, uses production function from Smith 2008; else,
uses hartmann 1 Ga density, times ratio, divided by age.

default inputs:
age = 1
CF_diamter_min = 0.0039
diameter_max = 1024
EROSION = false
SECONDARY = false
Beta = 1e-5
HARTMANN_PROD = true
"""
function craterfreq(age::Number=1, diameter_min::Number=0.0039,
    diameter_max::Number=1024, EROSION::Bool=false, SECONDARY::Bool=false,
    Beta::Number=1e-5; HARTMANN_PROD::Bool=true)

    ## Calculates (with erosion) or loads crater sizes and frequencies

    # Crater diameter bins, Michael Icarus 2013
    D = [0.00391, 0.00553, 0.00782, 0.0111, 0.01565, 0.0221, 0.0313, 0.0442,
            0.06251, 0.0887, 0.125, 0.177, 0.25, 0.354, 0.5, 0.7075, 1., 1.415,
            2., 2.83, 4., 5.66, 8.05, 11.32, 16.05, 22.63, 32.05, 45.3, 64.05,
            90.6, 128.05, 181.1, 256.05, 362.1, 512.05, 724.1,]  #km

    # Crater frequencies for 1 Ga surface, Michael Icarus 2013
    Nh = [4.04E+03, 2.33E+03, 1.14E+03, 4.58E+02, 1.91E+02, 6.66E+01, 2.40E+01,
             9.44E+00, 3.30E+00, 1.22E+00, 4.37E-01, 1.47E-01, 4.70E-02, 1.38E-02,
             4.02E-03, 1.15E-03, 3.08E-04, 1.28E-04, 6.85E-05, 3.67E-05, 1.98E-05,
             1.06E-05, 5.68E-06, 3.04E-06, 1.62E-06, 8.71E-07, 4.67E-07, 2.40E-07,
             1.12E-07, 5.21E-08, 2.43E-08, 1.13E-08, 5.28E-09, 2.47E-09, 1.15E-09,
             5.37E-10,]

    N_1Ga = 3.79e-14*(exp(6.93)-1)+5.84e-4   #Michael 2013 Icarus, pg 889, eqn 3
    ratio = (3.79e-14*(exp(6.93*age)-1)+5.84e-4*age)/N_1Ga  #ratio to 1 Ga

    if EROSION
        # Allocate array for results
        N = Array{Float64}(undef, size(Nh))

        # Set minimum erosion rate
        if Beta == 0 ## Q for Marisa: how about if Beta < 1e-5?
            Beta = 1e-5
        end

        for i = 1:length(D)
            # Initial conditions
            p = L = Ξ = Ψ = 0

            # Set production fnctn p
            if ~HARTMANN_PROD
               if D[i] < 1.4
                    p = 0.0035*(0.13*log(D[i])+0.83)/(D[i]^3.3)
                else
                    if D[i] <= 48.1
                        p = 10^(-1.8*log10(D[i])-2.59)
                    else
                        p = 10^(-2.2*log10(D[i])-1.89)
                    end
                end
                p = p * 0.29
            else
                p = Nh[i]*ratio/age
            end

            # Set Xi
            if D[i] < 5.8
                Ξ = 0.2*D[i]
            else
                Ξ = 0.42*log(D[i])-0.01
            end

            # Set Psi
            if D[i] < 1.2 && SECONDARY
                Ψ = 1 + 0.1225/(0.1225+D[i]^2)  # yields eqn 10
            else
                Ψ = 1
            end

            # lambda: crater loss fnctn
            L = Ψ*Beta/(1000*Ξ)

            N[i] = p/L*(1-exp(-L*age))
        end
    else
        N = Nh .* ratio
    end

    ## Crater vectors of crater sizes and frequencies based on desired crater size range

    # Find closest value in D to diameter_min
    CF_idx_min = argmin(abs.(D .- diameter_min))

    # Find closest value in D to diameter_max
    CF_idx_max = argmin(abs.(D .- diameter_max))

    # create arrays for diameter range and crater frequency range
    CF_crater_D = D[CF_idx_min:CF_idx_max]
    CF_crater_N = N[CF_idx_min:CF_idx_max]

    return (CF_crater_D, CF_crater_N)
end

## --- Define cumul function

function cumul(crater_count_density, IGNORE_EMPTY::Bool=true)
    ## Cumulates crater count densities
    num_bins = length(crater_count_density)
    cumulate_count_density = Array{Float64}(undef,num_bins)
    for k = 1:num_bins   #size bins
        cumulate_count_density[k] = 0
        if (crater_count_density[k] > 0) || ~IGNORE_EMPTY
            for m = 1:(num_bins - k + 1)
                cumulate_count_density[k] += crater_count_density[num_bins - m + 1]
            end
        end
    end

    return cumulate_count_density
end

## ---

function craterfit(diameters, count_dens, weight, CUMULATIVE::Bool,
    SKIP_EMPTY::Bool, FIT_EROSION::Bool=false, SECONDARY::Bool=false)

    ## Provides best fit age and beta (if FIT_EROSION) and RMSE for a crater count

    # Precision for ages: higher slows runtime but gives more precise ages
    FIT_prec = 10

    if FIT_EROSION  #Test betas if FIT_EROSION is on
        test_b = [1, (10:10:400)...]
    else
        test_b = [1]
    end

    if ~CUMULATIVE
        # Just the diameters measured
        max_d = diameters[end]
    else
        # All craters larger than measured, since need to cumulate them
        max_d = 1000
    end

    RMSE = Array{Float64,2}(undef, length(test_b), 4*FIT_prec)
    for i = 1:length(test_b)  # betas
        for k = 1:(4*FIT_prec) # ages

            RMSE[i,k] = 0

            if FIT_EROSION   #gets theoretical crater frequency
                d_theory = craterfreq(k/FIT_prec, diameters[1],
                    max_d, true, SECONDARY, test_b[i])[2]
            else
                d_theory = craterfreq(k/FIT_prec, diameters[1], max_d)[2]
            end

            if CUMULATIVE  #Cumulates crater frequnecy if desired
                d_theory = cumul(d_theory, SKIP_EMPTY)
            end

            for m = 1:length(diameters)
                if weight==4  #gets RMSE for log weighting
                    if count_dens[m]>0
                        RMSE[i,k] .= ( RMSE[i,k]
                            + (log10(count_dens[m])
                            - log10(d_theory[m]))^2 )
                    end
                else
                    if ~isnan(count_dens[m])  #gets RMSE for non-log weighting
                        if ~SKIP_EMPTY || (count_dens[m] != 0)
                            RMSE[i,k] = ( RMSE[i,k]
                                + ((count_dens[m] - d_theory[m]).^2
                                / d_theory[m].^(weight-1)) )
                        end
                    end
                end
            end
            RMSE[i,k] = RMSE[i,k]^0.5
        end
    end

    if FIT_EROSION  #finds minimum RMSE, depending on whether FIT_EROSION, which requires checking multiple columns
        RMSE_min = minimum(minimum(RMSE, dims=1))
        idx_k = argmin(minimum(RMSE, dims=1)[:])
        idx_i = argmin(RMSE[:,idx_k])
        age_calc = idx_k/FIT_prec
        beta_calc = test_b[idx_i]
    else
        RMSE_min = minimum(RMSE[:])
        idx = argmin(RMSE[:])
        age_calc = idx / FIT_prec
        beta_calc = 0
    end

    return (age_calc, beta_calc, RMSE_min)
end

## ---

# Plots calculated subsurface ages given a parent age, surface area, and Beta value
function CraterModel(;min_age::Number=1, max_age::Number=3, num_ages::Int=3,
    min_crat_matrix::Array{<:Number}=[.016, .063, .125, .178, .25, .375, .5, 1.],
    area_surface::Number=20000, sub_area=100, iters::Int=100, diameter_max::Number=1000,
    beta::Number=10, weight::Number=2, EROSION::Bool=true, SECONDARY::Bool=false,
    CUMULATIVE::Bool=true, SKIP_EMPTY::Bool=true, USE_QUANTILE::Bool=true,
    FIT_EROSION::Bool=true, AVG_NO_CRATERS::Bool=false)

    # # Defaults:
    # min_age = 1 # min age to be considered
    # max_age = 3  # max age to be considered
    # num_ages = 3 # number of ages to be considered
    # area_surface = 20000 # large area with known age (i.e., true_age_surface), in km2. 200000 used for subsurface areas greater than 10 km
    # sub_area = 100       # supsampled area, km
    # iters = 100  #number of iterations
    # diameter_max = 1000 # largest crater diameter on area, in km (max = 1024 km)
    # EROSION = true # if erosion is modeled Based on Smith
    # SECONDARY = false  # if secondary craters are modeled based on SMith
    # CUMULATIVE = true # counts and plots cumulatively
    # SKIP_EMPTY = true # skips bins with no craters when using cumulative
    # USE_QUANTILE = true  # u se interquartile range or mean/std dev to invert?
    # beta = 10 # erosion rate in nm/a
    # weight = 2
    # FIT_EROSION = true  # If fitting counted crater densities uses Smith model to include erosion
    # AVG_NO_CRATERS = false  # If should average subareas with no craters or ignore them (have 0.1 Ga age)

    # Create age and min crater size matrices
    age_matrix = collect(range(min_age, max_age, length=num_ages))

    # Alocate intermediate result arrays
    age_calc = Array{Float64}(undef, num_ages, length(sub_area), iters, length(min_crat_matrix))
    beta_calc = Array{Float64}(undef, num_ages, length(sub_area), iters, length(min_crat_matrix))
    RMSE_min = Array{Float64}(undef, num_ages, length(sub_area), iters, length(min_crat_matrix))

    # Allocate output arrays
    age_avg = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    age_std = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    beta_avg = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    beta_std = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    age_med = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    age_25 = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))
    age_75 = Array{Float64}(undef, num_ages, length(sub_area), length(min_crat_matrix))

    ## Calculates average ages
    @showprogress for h = 1:num_ages

        # crater density
        age = age_matrix[h]

        # find diameters and # of craters
        (diameters, density_craters) = craterfreq(age, min_crat_matrix[1],
            diameter_max, EROSION, SECONDARY, beta);

        area_crater = pi*(diameters/2).^2; #computes crater areas

        # Calculates total number of craters after rounding
        num_craters = round.(Int, area_surface .* density_craters)
        num_bins = length(diameters)
        total_craters = sum(num_craters)
        max_craters_per_bin = maximum(num_craters)

        # create matrix of (crater diameters bins, each crater in the bin, x&y)
        # create coordinates for edges of surface
        x_length = sqrt(area_surface)
        y_length = sqrt(area_surface)

        x_min = 0
        y_min = 0
        x_max = sqrt(area_surface)
        y_max = sqrt(area_surface)

        # address_book = Array{Float64}(undef, num_bins, max_craters_per_bin, 2)
        address_book = fill(NaN, num_bins, max_craters_per_bin, 2)

        for i=1:num_bins
            if num_craters[i] > 0   # checks if any craters are in this diameter bin
                address_book[i,1:num_craters[i],1] .= x_length .* rand(num_craters[i])
                address_book[i,1:num_craters[i],2] .= y_length .* rand(num_craters[i])
            end
        end

        # Counts craters located within sub-sample areas (each sub-sampled area
        # placed on the board randomly and number of sub-samples determined by
        # "iters"

        # counts number of craters for each sub area, for each iteration, in each bin
        crater_count = fill(0, length(sub_area), iters, num_bins)
        crater_count_density = Array{Float64}(undef, length(sub_area), iters, num_bins)
        box = Array{Float64}(undef, length(sub_area), iters, 2)

        for a = 1:length(sub_area)
            for b = 1:iters
                #Defines a random box for this iteration
                box[a,b,1] = rand() * (x_max - sqrt(sub_area[a]))
                box[a,b,2] = rand() * (y_max - sqrt(sub_area[a])) # coordinates of subsurface a

                #Checks each crater to see if it's in this interation's box

                for c = 1:num_bins # Each diameter bin size
                    # If there are craters in bin smaller than subsurface
                    if (num_craters[c] > 0) && (area_crater[c] < sub_area[a])
                        for d=1:num_craters[c] # Each crater in the bin
                            if (address_book[c,d,1] >= box[a,b,1]) && (address_book[c,d,1] <= (box[a,b,1]+sqrt(sub_area[a]))) &&
                                (address_book[c,d,2] >= box[a,b,2]) && (address_book[c,d,2] <= (box[a,b,2]+sqrt(sub_area[a])))

                                crater_count[a,b,c] += 1
                            end
                        end
                    end
                    crater_count_density[a,b,c] = crater_count[a,b,c]/sub_area[a]
                end
            end
        end

        # Calculates RMSE compared to isochrons with no erosion
        for k = 1:length(min_crat_matrix)
            idx = argmin(abs.(diameters .- min_crat_matrix[k])) # bin for min crater diameter
            for i = 1:length(sub_area) # subsurface areas
                for j = 1:iters  # iterations
                    if CUMULATIVE # cumulates crater densities if desired
                        RMSE_density = cumul(crater_count_density[i,j,:], SKIP_EMPTY)
                    else
                        RMSE_density = crater_count_density[i,j,:]
                    end

                    (age_calc[h,i,j,k], beta_calc[h,i,j,k], RMSE_min[h,i,j,k]) =
                        craterfit(diameters[idx:end], RMSE_density[idx:end], weight, CUMULATIVE, SKIP_EMPTY, FIT_EROSION, SECONDARY)
                end
                # Calculates averages and standard deviations over all iterations
                if AVG_NO_CRATERS
                    age_to_avg = age_calc[h,i,:,k]  #includes all ages, even if no craters (0.1 Ga)
                    beta_to_avg = beta_calc[h,i,:,k]
                else
                    age_to_avg = age_calc[h,i,age_calc[h,i,:,k].>0.1,k]  #only includes ages with craters (>0.1 Ga)
                    beta_to_avg = beta_calc[h,i,beta_calc[h,i,:,k].>0.1,k]
                end
                age_avg[h,i,k] = isempty(age_to_avg) ? NaN : mean(age_to_avg)  #calculates mean and std dev to plot
                age_std[h,i,k] = isempty(age_to_avg) ? NaN : std(age_to_avg)
                beta_avg[h,i,k] = isempty(beta_to_avg) ? NaN : mean(beta_to_avg)
                beta_std[h,i,k] = isempty(beta_to_avg) ? NaN :  std(beta_to_avg)
                age_med[h,i,k] = isempty(age_to_avg) ? NaN : median(age_to_avg)  #calculates median and interquantile to plot
                age_25[h,i,k] = isempty(age_to_avg) ? NaN : percentile(age_to_avg, 25)
                age_75[h,i,k] = isempty(age_to_avg) ? NaN : percentile(age_to_avg, 75)
            end
        end
    end
    return (age_matrix, age_avg, age_std, beta_avg, beta_std, age_med, age_25, age_75)
end

## --- Run and plot results

    min_crat_matrix=[.016, .063, .125, .178, .25, .375, .5, 1.]
    sub_area = [100]
    (age_matrix, age_avg, age_std, beta_avg, beta_std, age_med, age_25, age_75) = CraterModel(min_crat_matrix=min_crat_matrix, sub_area=sub_area)


    f = plot(xlabel="Minimum crater diameter (m)", ylabel="Median calculated age (Ga)",
        fg_color_legend=:white, legend=:bottomright)
    for i = 1:length(age_matrix)
        med = age_med[i,1,:]
        quant25 = med .- age_25[i,1,:]
        quant75 = age_75[i,1,:] .- med
        plot!(f, min_crat_matrix.*1000, med, yerror=(quant25, quant75),
            seriestype=:scatter, markershape=:auto, mscolor=:auto,
            label="$(age_matrix[i])")
    end
    plot!(f, title="$(sub_area[1]) km^2", xscale=:log10, ylims=(0,4), xlims=(10,1200))
    plot!(f, xticks=([16,63,125,250,500,1000],["16","63","125","250","500","1000"]))
    display(f)


### ---

## --- Generate Figure 2 or 5

# Inputs
min_crat_matrix=[.016, .063, .125, .178, .25, .375, .5, 1.]
sub_area = 100

# Seed the global rng
Random.seed!(reinterpret(Int,time()))

# Run the model
(age_matrix, age_avg, age_std, beta_avg, beta_std, age_med, age_25, age_75) =
    CratModelFig2and5(min_crat_matrix=min_crat_matrix, sub_area=sub_area)

f = plot(xlabel="Minimum crater diameter (m)", ylabel="Median calculated age (Ga)",
    fg_color_legend=:white, legend=:bottomright)
for i = 1:length(age_matrix)
    med = age_med[i,:]
    quant25 = med .- age_25[i,:]
    quant75 = age_75[i,:] .- med
    plot!(f, min_crat_matrix.*1000, med, yerror=(quant25, quant75),
        seriestype=:scatter, markershape=:auto, mscolor=:auto,
        label="$(age_matrix[i])")
end
plot!(f, title="$(sub_area) km^2", xscale=:log10, ylims=(0,4), xlims=(10,1200))
plot!(f, xticks=([16,63,125,250,500,1000],["16","63","125","250","500","1000"]))
savefig(f, "age_diameter_scaling.pdf")
display(f)

## --- Generate Figure 7 or

filepath = "QuantainLandslideCraters.diam"
CUMULATIVE = false
SKIP_EMPTY = true
iters = 100
prec = 0.25

# Seed the global rng
Random.seed!(reinterpret(Int,time()))

(fit_ages, fit_betas, fit_craters, diameters, compare_count, sub_surf) =
    CratModelFig7and8(filepath=filepath, CUMULATIVE=CUMULATIVE, SKIP_EMPTY=SKIP_EMPTY,
                      iters=iters, prec=prec)

age_hist_bins = 0:0.5:4
beta_hist_bins = 0:50:400
diameter_min = minimum(diameters)
plottitle = "$filepath; $sub_surf km^2; Dmin: $diameter_min km;$(CUMULATIVE ? "" : " not") cumulative; iters: $iters; prec: $prec"

## ---  2d-histogram (age and beta)

h1 = plot(xlabel="Age [Ga]", ylabel="\\beta [{nm/a}]", xlims=(0,4), ylims=(0,400), title=plottitle, titlefontsize=8)
histogram2d!(h1, fit_ages, fit_betas, fc=:plasma, bins=(age_hist_bins,beta_hist_bins), colorbar_title="Modeled surfaces")
savefig(h1, "craterhistogram_age_beta.png")
display(h1)

## --- 1-d histogram, marginalized on age

w = fit(Histogram,fit_ages,age_hist_bins).weights
binwidth = abs(age_hist_bins[2]-age_hist_bins[1])
h2 = plot(xlabel="Age [Ga]", ylabel="Percent of matching surfaces", xlims=(0,4), ylims=(0,100))
bar!(h2, cntr(age_hist_bins), w*100/sum(w), label="", bar_width=binwidth, linealpha=0)
plot!(h2, yticks=(0:10:100, (0:10:100).|> x-> "$x %"), title=plottitle, titlefontsize=8)
savefig(h2, "craterhistogram_marginalized_age.pdf")
display(h2)
## --- Compare modelled and measured counts

# Which fitted count to display
plot_num = 1
plot_data = fit_craters[plot_num,:]./sub_surf |> x -> CUMULATIVE ? cumul(x,SKIP_EMPTY) : x

# Make the plot
h3 = plot(xscale=:log10, yscale=:log10, xlabel="Crater Diameter [km]", fg_color_legend=:white)
plot!(h3, xlims=(10^-3,10^3), ylims=(10^-6, 10^3))
plot!(h3, diameters, plot_data, label="Modeled", seriestype=:scatter)
plot!(h3, diameters, compare_count./sub_surf, label="Measured", seriestype=:scatter, shape=:star5)
plot!(h3, ylabel = "C$(CUMULATIVE ? "umulative c" : "")rater count density ($(CUMULATIVE ? "cumulative " : "")craters / km^2)")
plot!(h3, title="$(fit_ages[plot_num]) Ga; \\beta = $(fit_betas[plot_num]) nm a^{-1}; $plot_num of $(length(fit_ages))" )
savefig(h3, "model-data_comparison_$(plot_num).pdf")
display(h3)

## ----
