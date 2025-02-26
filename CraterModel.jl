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

    function crater_mcmc(diam, density, density_sigma; 
            nsteps=100000, 
            burnin=nsteps, 
            age=3.5, 
            erosion=1.0,
            age_prior = Uniform(0, 4.567),
            erosion_logprior = Normal(0,4),
        )

        # Allocate variables to record stationary distribution
        acceptancedist = falses(nsteps)
        lldist = zeros(nsteps)
        agedist = zeros(nsteps)
        erosiondist = zeros(nsteps)

        σj_agedist = zeros(nsteps)
        σj_erosiondist = zeros(nsteps)

        # Initial proposal
        ageₚ = age
        erosionₚ = erosion
        densityₚ = copy(density)
        llₚ = ll = crater_ll!(densityₚ, diam, density, density_sigma, age, erosion) + logpdf(age_prior, age) + logpdf(erosion_logprior, log(erosion))

        # Jumping distributions
        σj_age = 0.04
        σj_logerosion = log(2)

        # Run Markov chain!
        for i in 1:nsteps+burnin
            # Make new proposal based on last accepted propsoal
            r = rand(1:2)
            if r == 1 
                ageₚ = age + rand(Normal(0, σj_age))
            elseif r == 2
                erosionₚ = exp(log(erosion) + rand(Normal(0,σj_logerosion)))
            end

            # Calculate log likelihood of new proposal
            llₚ = crater_ll!(densityₚ, diam, density, density_sigma, ageₚ, erosionₚ) + logpdf(age_prior, ageₚ) + logpdf(erosion_logprior, log(erosionₚ))
            
            # Accept or reject proposal
            if log(rand()) < (llₚ-ll)
                # Update jumping distributions
                ageₚ ≠ age && (σj_age = 2.5 * abs(ageₚ - age))
                erosionₚ ≠ erosion && (σj_logerosion = 2.5 * abs(log(erosionₚ)-log(erosion)))

                # Update accepted values
                ll = llₚ
                age = ageₚ
                erosion = erosionₚ
                if i > burnin
                    acceptancedist[i-burnin] = true
                end
            end

            # Record current accepted result
            if i > burnin
                lldist[i-burnin] = ll
                agedist[i-burnin] = age
                erosiondist[i-burnin] = erosion

                σj_agedist[i-burnin] = σj_age
                σj_erosiondist[i-burnin] = σj_logerosion
            end
        end
        @info "mean acceptance: $(mean(acceptancedist))\n"*
              "  σj_age: $(round(mean(σj_agedist), sigdigits=6))\n"*
              "  σj_logerosion: $(round(mean(σj_erosiondist), sigdigits=6))\n"

        return acceptancedist, lldist, agedist, erosiondist
    end

    function crater_ll!(densityₚ, diam, density, density_sigma, age, erosion)
        diam_cf, density_cf = craterfreq(age, erosion)
        # Interpolate, in log space
        density_cf .= max.(density_cf, 0.0)
        density_cf .= log.(density_cf)
        linterp1!(densityₚ, diam_cf, density_cf, diam)
        densityₚ .= exp.(densityₚ)
        # Calculate log likelihood
        ll = 0.0
        for i in eachindex(densityₚ, density, density_sigma)
            if !isnan(density[i])
                ll += logpdf(Normal(density[i], density_sigma[i]), densityₚ[i]) # Could we use Poisson directly here??
                # logpdf(Normal(mean, sigma), comparison?) = log likelihood of observing densityp assuming density is normally distributed around density
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


