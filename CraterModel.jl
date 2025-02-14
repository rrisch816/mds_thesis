## --- Load required packages

    using Statistics, StatsBase, Random, LaTeXStrings, Plots
    using ProgressMeter: @showprogress
    if ~ @isdefined cntr # if cntr not defined,
        function cntr(edges) # computes bin center from edges
            return (edges[1:end-1]+edges[2:end])/2
        end
    end
    function nextnicenumber(x::Number) # finds next nice number greater or equal to x
        nicenumber = [1,1.2,1.25,1.4,1.5,1.6,1.75,1.8,2,2.5,3,4,5,6,7.5,8,9,10,]
        # nicedivisor = [5, 6,   5,  7,  5,  4,   7,  6,4,  5,3,4,5,6,  5,4,3, 5,]
        base = 10^floor(log10(x)) # takes log base 10 of x (power of 10 closest to x), extracts integer part of the exponent, reconstructs nearest power of 10 <= x
        i = findfirst(nicenumber .>= x/base) # scale x down by dividing by base, returns index of that number
        return nicenumber[i]*base
    end

## --- Add MCMC inversion function!

    function crater_mcmc(diam, density, density_sigma; nsteps=100000, burnin=1000, age=3.5, erosion=1e1)
        # Allocate variables to record stationary distribution
        lldist = zeros(nsteps)
        agedist = zeros(nsteps)
        erosiondist = zeros(nsteps)

        # Initial proposal
        ageₚ = age
        erosionₚ = erosion
        densityₚ = copy(density)
        llₚ = ll = crater_ll!(densityₚ, diam, density, density_sigma, age, erosion)

        for i in 1:nsteps+burnin
            # Make new proposal based on last accepted propsoal
            erosionₚ = exp(log(erosion) + rand(Normal(0,1)))
            ageₚ = age + rand(Normal(0, 0.05))

            # Calculate log likelihood of new proposal
            llₚ = crater_ll!(densityₚ, diam, density, density_sigma, ageₚ, erosionₚ)
            
            if log(rand()) < (llₚ-ll)
                # Accept proposal!
                ll = llₚ
                age = ageₚ
                erosion = erosionₚ
            end

            # Record current accepted result
            if i > burnin
                lldist[i-burnin] = ll
                agedist[i-burnin] = age
                erosiondist[i-burnin] = erosion
            end
        end

        return lldist, agedist, erosiondist
    end

    function crater_ll!(densityₚ, diam, density, density_sigma, age, erosion)
        diam_cf, density_cf = craterfreq(age, erosion)
        # Possibly switch to interpolating in log space?
        linterp1!(densityₚ, diam_cf, density_cf, diam)
        ll = 0.0
        for i in eachindex(densityₚ, density, density_sigma)
            if !isnan(density[i])
                ll += logpdf(Normal(density[i], density_sigma[i]), densityₚ[i]) # Could we use Poisson directly here??
            end
        end
        return ll
    end

## --- Define craterfreq function

"""
'''julia
function craterfreq(age::Number=1, beta::Number=1e-5;
    \tdiameter_min::Number = 0.0039,
    \tdiameter_max::Number = 1024, 
    \tEROSION::Bool = beta>0, 
    \tSECONDARY::Bool=false,
    \tHARTMANN_PROD::Bool=true,
)
```

Based on range of craters chosen and surface age (Ga), returns list of crater
frequencies for 1 Ga surface age (no. craters/km2). If EROSION is true uses
Smith et al 2008 model for erosion and B is erosion rate (nm/a).

If HARTMANN_PROD if false, uses production function from Smith 2008; else,
uses hartmann 1 Ga density, times ratio, divided by age.
"""
function craterfreq(age::Number=1, beta::Number=1e-5;
        diameter_min::Number = 0.0039,
        diameter_max::Number = 1024, 
        EROSION::Bool = beta>0, 
        SECONDARY::Bool=false,
        HARTMANN_PROD::Bool=true,
    )
    # age defaults to 1 Ga
    # diameter_min and diameter_max defines crater diameter range
    # erosion if true, erosion effects applied (Smith et al 2008 model for erosion, beta is erosion rate, min 1e-5)
    # if secondary impact is true, use the secondary correction model
    # hartmann prod if true, uses hartman production function

    ## Calculates (with erosion) or loads crater sizes and frequencies

    # Crater diameter bins, Michael Icarus 2013
    # predefined logarithmically to cover wide range of diameters
    D = [0.00391, 0.00553, 0.00782, 0.0111, 0.01565, 0.0221, 0.0313, 0.0442,
            0.06251, 0.0887, 0.125, 0.177, 0.25, 0.354, 0.5, 0.7075, 1., 1.415,
            2., 2.83, 4., 5.66, 8.05, 11.32, 16.05, 22.63, 32.05, 45.3, 64.05,
            90.6, 128.05, 181.1, 256.05, 362.1, 512.05, 724.1,]  #km

    # Crater frequencies for 1 Ga surface, Michael Icarus 2013
    # Nh[i] is number of craters per km^2 for each diameter bin D[i]
    Nh = [4.04E+03, 2.33E+03, 1.14E+03, 4.58E+02, 1.91E+02, 6.66E+01, 2.40E+01,
             9.44E+00, 3.30E+00, 1.22E+00, 4.37E-01, 1.47E-01, 4.70E-02, 1.38E-02,
             4.02E-03, 1.15E-03, 3.08E-04, 1.28E-04, 6.85E-05, 3.67E-05, 1.98E-05,
             1.06E-05, 5.68E-06, 3.04E-06, 1.62E-06, 8.71E-07, 4.67E-07, 2.40E-07,
             1.12E-07, 5.21E-08, 2.43E-08, 1.13E-08, 5.28E-09, 2.47E-09, 1.15E-09,
             5.37E-10,]

    # Multiplicative factor to scale chronology function
    N_1Ga = 3.79e-14*(exp(6.93)-1)+5.84e-4   #Michael 2013 Icarus, pg 889, eqn 3

    # adjusts 1 Ga crater frequencies to match chosen surface age (age)
    ratio = (3.79e-14*(exp(6.93*age)-1)+5.84e-4*age)/N_1Ga  #ratio to 1 Ga

    # if Erosion is true
    if EROSION
        # empty array for adjusted crater frequencies
        N = Array{Float64}(undef, size(Nh))

        # Set minimum erosion rate
        if beta == 0 ## Q for Marisa: how about if beta < 1e-5?
            beta = 1e-5
        end

        #iterates through each crater diameter BIN to calculate crater loss due to erosion
        for i = 1:length(D)
            # Initial conditions
            p = L = Ξ = Ψ = 0

            # if not using hartmanns model, apply alt log production fcn
            # Set production fnctn p
            if ~HARTMANN_PROD
               if D[i] < 1.4 # for small craters
                    p = 0.0035*(0.13*log(D[i])+0.83)/(D[i]^3.3)
                else # for medium craters 
                    if D[i] <= 48.1
                        p = 10^(-1.8*log10(D[i])-2.59)
                    else # large craters
                        p = 10^(-2.2*log10(D[i])-1.89)
                    end
                end
                p = p * 0.29
            else # otherwise, if hartmanns is true, scales Nh[i] using age ratio
                p = Nh[i]*ratio/age
            end

            # Compute Depth Correction
            # Set Xi
            if D[i] < 5.8
                Ξ = 0.2*D[i] # correction factor, accounts for depth-dependent crater erosion
            else # large craters erode differently than small ones
                Ξ = 0.42*log(D[i])-0.01
            end

            # Secondary crater correction (ejecta impacts)
            # Set Psi
            if D[i] < 1.2 && SECONDARY
                Ψ = 1 + 0.1225/(0.1225+D[i]^2)  # yields eqn 10
            else
                Ψ = 1
            end

            # lambda: crater loss fnctn (rate of crater removal due to erosion (erosion rate beta))
            L = Ψ*beta/(1000*Ξ) 

            # adjusted crater frequency
            N[i] = p/L*(1-exp(-L*age))
        end
    else # if erosion false, simple scaling by ratio (Revisit, is redundant?)
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

    # gives x bins (diameter), y bins (frequencies)
    return (CF_crater_D, CF_crater_N)
end

## End RR commenting for now
# experiment with crater stats (does not account for erosion, rolloff)

# try plotting a couple units- those frequencies vs crater freq
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

## --- Define craterfit function

# repalce with MCMC?
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

## --- Define CratModelFig2and5 function

# Plots calculated subsurface ages given a parent age, surface area, and beta value
function CratModelFig2and5(;min_age::Number=1, max_age::Number=3, num_ages::Int=3,
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
    age_calc = Array{Float64}(undef, num_ages, iters, length(min_crat_matrix))
    beta_calc = Array{Float64}(undef, num_ages, iters, length(min_crat_matrix))
    RMSE_min = Array{Float64}(undef, num_ages, iters, length(min_crat_matrix))

    # Allocate output arrays
    age_avg = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    age_std = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    beta_avg = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    beta_std = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    age_med = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    age_25 = Array{Float64}(undef, num_ages, length(min_crat_matrix))
    age_75 = Array{Float64}(undef, num_ages, length(min_crat_matrix))


    sub_area_width = sqrt(sub_area)

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


        # Count craters located within sub-sample area placed on the board
        # randomly, with number of sub-samples determined by "iters"

        crater_count = fill(0, iters, num_bins)
        crater_count_density = Array{Float64}(undef, iters, num_bins)

        # count number of craters in each bin for each iteration
        for c = 1:num_bins # Each diameter bin size
            crater_x = x_length .* rand(num_craters[c])
            crater_y = y_length .* rand(num_craters[c])

            # If there are craters in bin smaller than subsurface
            if (num_craters[c] > 0) && (area_crater[c] < sub_area)
                for b = 1:iters
                    # Defines a random box for this iteration
                    box_x = rand() * (x_max - sub_area_width) # x-coordinate of subsurface area
                    box_y = rand() * (y_max - sub_area_width) # y-coordinate of subsurface area

                    # Checks each crater to see if it's in this iteration's box
                    @inbounds @simd for d=1:num_craters[c] # Each crater in the bin
                        if (crater_x[d] > box_x) && (crater_x[d] < (box_x+sub_area_width)) &&
                            (crater_y[d] > box_y) && (crater_y[d] < (box_y+sub_area_width))

                            crater_count[b,c] += 1
                        end
                    end
                end
            end
        end
        crater_count_density = crater_count ./ sub_area

        # Calculates RMSE compared to isochrons with no erosion
        for k = 1:length(min_crat_matrix)
            idx = argmin(abs.(diameters .- min_crat_matrix[k])) # bin for min crater diameter
            for j = 1:iters  # iterations
                if CUMULATIVE # cumulates crater densities if desired
                    RMSE_density = cumul(crater_count_density[j,:], SKIP_EMPTY)
                else
                    RMSE_density = crater_count_density[j,:]
                end

                (age_calc[h,j,k], beta_calc[h,j,k], RMSE_min[h,j,k]) =
                    craterfit(diameters[idx:end], RMSE_density[idx:end], weight, CUMULATIVE, SKIP_EMPTY, FIT_EROSION, SECONDARY)
            end
            # Calculates averages and standard deviations over all iterations
            if AVG_NO_CRATERS
                age_to_avg = age_calc[h,:,k]  #includes all ages, even if no craters (0.1 Ga)
                beta_to_avg = beta_calc[h,:,k]
            else
                age_to_avg = age_calc[h,age_calc[h,:,k].>0.1,k]  #only includes ages with craters (>0.1 Ga)
                beta_to_avg = beta_calc[h,beta_calc[h,:,k].>0.1,k]
            end
            age_avg[h,k] = isempty(age_to_avg) ? NaN : mean(age_to_avg)  #calculates mean and std dev to plot
            age_std[h,k] = isempty(age_to_avg) ? NaN : std(age_to_avg)
            beta_avg[h,k] = isempty(beta_to_avg) ? NaN : mean(beta_to_avg)
            beta_std[h,k] = isempty(beta_to_avg) ? NaN :  std(beta_to_avg)
            age_med[h,k] = isempty(age_to_avg) ? NaN : median(age_to_avg)  #calculates median and interquantile to plot
            age_25[h,k] = isempty(age_to_avg) ? NaN : percentile(age_to_avg, 25)
            age_75[h,k] = isempty(age_to_avg) ? NaN : percentile(age_to_avg, 75)
        end
    end
    return (age_matrix, age_avg, age_std, beta_avg, beta_std, age_med, age_25, age_75)
end

## --- Define CratModelFig7and8

#Reads .diam file and compares it to modeled surfaces from different ages
function CratModelFig7and8(;filepath::String="QuantainLandslideCraters.diam",
    ages::Array{<:Number}=collect(0.1:0.1:4.0), betas::Array{<:Number}=collect(0.0:8.0:400.0),
    prec::Number=0.25, iters::Int=100, area_surface::Number=10000,
    EROSION::Bool=true, SECONDARY::Bool=true, CUMULATIVE::Bool=false,
    SKIP_EMPTY::Bool=true, USE_COUNT_MIN_MAX::Bool=true,
    diameter_min::Number=0.088, diameter_max::Number=1000)

    # # Defaults:
    # ages = collect(0.1:0.1:4.0)    # Age bins
    # betas = collect(0.0:8.0:400.0) # Erosion rates
    # prec = 0.25;          #precision to consider crater count densities the same. 0.25 seems good.
    # iters = 100;          #iterations
    # area_surface = 10000; # large area with known age (i.e., true_age_surface), in km2
    # SECONDARY = true;  #apply Smith's secondary crater model or not--true for fig 7 false for fig 8
    # EROSION = true;     #whether or not to use erosion model for surface. should be on
    # CUMULATIVE = false; #cumulative crater plots--false for fig 7 true for fig 8
    # SKIP_EMPTY = true;  #skip empty bins with cumulative
    # USE_COUNT_MIN_MAX = true;  #whether or not to use the max/min bins from input .diam file
    #                             #true for fig 7 false for fig 8
    # diameter_min = .088;      #used if USE_COUNT_MIN_MAX is false
    # diameter_max = 1000;    #used if USE_COUNT_MIN_MAX is false

    # Length of input arrays
    num_ages = length(ages)
    num_betas = length(betas)

    # Open diam file and reads craters
    l = "start"
    craters = Array{Float64}(undef,0)
    surf_age = surf_area = surf_beta = sub_surf = nothing
    fid = open(filepath)
    while ~isempty(l)
        l = readline(fid)
        if occursin("age", l)
            surf_age = tryparse(Float64, match(r":(.*)$", l)[1])
        elseif occursin("Large", l)
            surf_area = tryparse(Float64, match(r":(.*)$", l)[1])
        elseif occursin("B", l)
            surf_beta = tryparse(Float64, match(r"=sur(.*)$", l)[1])
        elseif occursin("area", lowercase(l))
            sub_surf = tryparse(Float64, match(r"=(.*)$", l)[1])
        elseif occursin("crater = {diameter", l)
            p = tryparse(Float64, readline(fid))
            while ~isnothing(p)
                push!(craters, p)
                p = tryparse(Float64, readline(fid))
            end
        end
    end
    close(fid)

    # Set diameter bins, from 0.0039-724 km
    diameters = craterfreq(1,0,1200)[1]

    # Bin craters
    compare_crater_count = Array{Int64}(undef, length(diameters))
    for i = 1:(length(diameters)-1)
        compare_crater_count[i] = count((craters .>= diameters[i]) .& (craters .< diameters[i+1]))
    end
    compare_crater_count[end] = count(craters .>= diameters[end])  # last bin

    if USE_COUNT_MIN_MAX
        idx_min = findfirst(compare_crater_count .> 0)      #finds smallest bin
        idx_max = findlast(compare_crater_count .> 0)      #finds largeest bin
        diameters = diameters[idx_min:idx_max]        #takes smallest to largest present bins
        compare_crater_count = compare_crater_count[idx_min:idx_max]  #takes smallest to largest present bins
    else  # user-defined crater sizes
        idx_min = argmin(abs.(diameter_min .- diameters))
        idx_max = argmin(abs(diameter_max .- diameters))
        diameters = diameters[idx_min:idx_max]
        compare_crater_count = compare_crater_count[idx_min:idx_max]
    end

    if CUMULATIVE  #cumulates crates if needed
        compare_crater_count .= cumul(compare_crater_count, SKIP_EMPTY)
    end

    # Set diameter min and max to exact values from filtered array
    diameter_min = minimum(diameters)
    diameter_max = maximum(diameters)
    num_bins = length(diameters)

    # Assign output variables -- maximum possible size
    fit_ages = Array{Float64}(undef, 0)
    fit_betas = Array{Float64}(undef, 0)
    fit_craters = Array{Float64}(undef, 0)
    test_crater_count = Array{Int64}(undef, num_bins)
    # perc_diff = Array{Float64}(undef, num_ages, num_betas, iters, num_bins)

    # Create and subsamples surfaces, and compares to inputted crater count
    sub_area_width = sqrt(sub_surf)
    fit_count = 0
    @showprogress for h = 1:num_ages
        age = ages[h]
        for k = 1:num_betas
            beta = betas[k] # Erosion rate to model

            # Find diameters and # of craters
            (test_diameters, test_density_craters) =
                craterfreq(age,diameter_min,diameter_max,EROSION,SECONDARY,beta)

            area_crater = pi * (test_diameters / 2).^2 #computes crater areas

            #Calculates total number of craters after rounding
            num_craters = round.(Int, area_surface.*test_density_craters)
            # num_bins = length(test_diameters)
            total_craters = sum(num_craters)
            max_craters_per_bin = maximum(num_craters)

            # Create coordinates for edges of surface
            x_length = sqrt(area_surface)
            y_length = sqrt(area_surface)
            x_max = sqrt(area_surface)
            y_max = sqrt(area_surface)

            # Random crater positions
            address_book_x = num_craters .|> x -> x_length .* rand(x)
            address_book_y = num_craters .|> x -> y_length .* rand(x)

            # Count number of craters for each sub area, for each iteration, in each bin
            for b=1:iters
                test_crater_count .= 0
                #Defines a random box for this iteration
                box_x = rand()*(x_max - sub_area_width)
                box_y = rand()*(y_max - sub_area_width)

                #Checks each crater to see if it's in this interation's box
                for c = 1:num_bins      #each diameter bin size
                    if (num_craters[c] > 0) && (area_crater[c] < sub_surf)    #If there are craters in bin smaller than subsurface
                        @inbounds @simd for d=1:num_craters[c]     #for each crater in the bin
                            if (address_book_x[c][d] > box_x) && (address_book_x[c][d] < (box_x+sub_area_width)) &&
                                (address_book_y[c][d] > box_y) && (address_book_y[c][d] < (box_y+sub_area_width))

                                test_crater_count[c] += 1
                            end
                        end
                    end
                end

                if CUMULATIVE  #cumulates synthetic count if needed
                    test_crater_count .= cumul(test_crater_count, SKIP_EMPTY)
                end

                # calculates difference between .diam count and synthetic count
                rel_diff = 2 * abs.(test_crater_count - compare_crater_count) ./
                    (test_crater_count + compare_crater_count) .|>
                    x -> isnan(x) ? 0 : x

                # perc_diff[h,k,b,:] .= rel_diff
                # if all(perc_diff[h,k,b,:] .< prec)  #checks that all bins meet prec criteria
                if all(rel_diff .< prec)  #checks that all bins meet prec criteria
                    fit_count += 1  # Counts number of fits
                    # fit_ages[fit_count] = age
                    # fit_betas[fit_count] = beta
                    # fit_craters[fit_count,:] .= test_crater_count
                    push!(fit_ages, age)   # Store age of fit
                    push!(fit_betas, beta) # Store beta of fit
                    append!(fit_craters, test_crater_count) # Store synthetic count for displaying later
                end
            end
        end
    end
    return (fit_ages, fit_betas, collect(reshape(fit_craters, 8, fit_count)'), diameters, compare_crater_count, sub_surf)
end

