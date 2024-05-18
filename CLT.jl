using Distributions, Random, DataFrames, CSV, StatsBase


# TODO replace with zscore from StatsBase
function standardize(x, μ, σ, n::Int64)::Float64
    (x - μ) / (σ / sqrt(n))
end

function t_score(data::Vector{Float64}; μ::Real)::Float64
    standardize(mean(data), μ, std(data), length(data))
end

function sampling_distribution(d::Distribution, n::Int, r::Int, statistic::Function; args...)::Vector{Float64}
    Random.seed!(0)
    # Preallocating vectors for speed
    sample_statistics = zeros(r) # TODO move out of this function
    sample = zeros(n)

    # samples = rand(d, n, r)
    # sample_statistics = mapcols(col -> statistic(col; args...), samples)

    # sample = rand(d, r, n) # Generate the sample matrix directly

    # # Define a helper function to apply statistic with optional arguments to each column
    # function apply_statistic(column)
    #     if isempty(args)
    #         return statistic(column)
    #     else
    #         return statistic(column, args...)
    #     end
    # end

    # # Apply the helper function to every column in the sample matrix
    # result = mapcols(apply_statistic, sample)

    # sample = zeros(r, n)
    # sample_statistics = statistic.(rand!(d, sample); args...)

    # Sampling r times and calculating the statistic
    @inbounds for i in 1:r
        rand!(d, sample) # in-place to reduce memory allocation
        sample_statistics[i] = statistic(sample; args...)::Float64

        # # TODO Hacky, fix
        # try
        #     sample_statistics[i] = statistic(sample; args...)::Float64
        # catch e
        #     if isa(e, MethodError)
        #         sample_statistics[i] = statistic(sample)::Float64
        #     else
        #         rethrow(e)
        #     end
        # end
    end


    #print dimensions of sample_statistics
    # println(size(sample_statistics))
    sample_statistics
end

function analysis(statistic::Function, d::Distribution, n::Int, r::Int, μ::Real, σ::Real, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64}
    sample_statistics = sampling_distribution(d, n, r, statistic; args...)
    skewness = StatsBase.skewness(sample_statistics)
    kurtosis = StatsBase.kurtosis(sample_statistics)

    # Standardizing the values to look at tail probabilities
    z_scores = standardize.(sample_statistics, μ, σ, n)

    # TODO Fix critical value
    upper = sum(z_scores .>= critical) / r
    lower = sum(z_scores .<= -critical) / r

    (upper, lower, skewness, kurtosis)
end

function analyze_distributions(statistic::Function, r::Int64, sample_sizes::Vector{Int64}, critical::Function, distributions, params=false)::DataFrame
    println("Analyzing sampling distributions of $(statistic) with $(r) repetitions")
    # Setting up the results we're interested in
    results = DataFrame(
        "Distribution" => String[],
        "Skewness" => Float64[],
        "Sample Size" => Int64[],
        "Upper Tail" => Float64[],
        "Lower Tail" => Float64[],
        "Sampling Skewness" => Float64[],
        "Sampling Kurtosis" => Float64[],
        "Population Mean" => Float64[],
        "Population SD" => Float64[]
    )

    # Analyzing each distribution
    @inbounds for d::Distribution in distributions
        println(string(d))

        # Getting population parameters
        μ = mean(d)
        σ = std(d)
        skewness = StatsBase.skewness(d)

        # Analyzing each sample size
        u = Threads.SpinLock() # lock to avoid data races
        @inbounds Threads.@threads for n in sample_sizes
            if params
                #TODO fix this
                upper, lower, sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)), μ=μ)
            else
                upper, lower, sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)))
            end
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, upper, lower, sample_skewness, sample_kurtosis, μ, σ))
            end
        end

    end

    sort!(results, [:Distribution, :"Sample Size"])
end

function main(r)
    sample_sizes = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500]
    distributions = [
        Gamma(16),
        LogNormal(0, 0.25),
        Gamma(4),
        Gamma(2),
        LogNormal(0, 0.5),
        Gamma(1),
        Exponential(),
        Gamma(0.64),
        LogNormal(0, 0.75)
    ]
    tstar = n -> quantile(TDist(n), 0.975)
    zstar = n -> quantile(Normal(), 0.975)

    # TODO add back in
    # Compile
    analyze_distributions(mean, 1, sample_sizes, zstar, distributions)
    analyze_distributions(t_score, 1, sample_sizes, tstar, distributions, true)

    # Warning: this code will take a very long time to run if used with a large r. We used r = 10_000_000
    @time means::DataFrame = analyze_distributions(mean, r, sample_sizes, tstar, distributions)
    CSV.write("means.csv", means)

    @time t::DataFrame = analyze_distributions(t_score, r, sample_sizes, tstar, distributions, true)
    CSV.write("t.csv", t)
end

function graphing(r)
    d = Exponential()
    μ = mean(d)
    σ = std(d)
    n = 30

    # Creating sampling distributions
    exponential30 = sampling_distribution(mean, d, n, r)
    exponential30std = standardize.(exponential30, μ, σ, n)

    n = 150
    exponential150 = sampling_distribution(mean, d, n, r)
    exponential150std = standardize.(exponential150, μ, σ, n)

    graphing = DataFrame(
        "Exponential" => exponential30,
        "Exponential Z-Scores" => exponential30std,
        "Exponential 150" => exponential150,
        "Exponential 150 Z-Scores" => exponential150std
    )

    CSV.write("graphing.csv", graphing)
end

# @profview main(100_000)
main(10_000_000)
# graphing(10_000_000)

#=
CORRECT T* VALUES
5 2.570581835636314
10 2.228138851986274
20 2.0859634472658644
30 2.0422724563012378
40 2.0210753903062737
50 2.0085591121007615
60 2.00029782201426
70 1.9944371117711865
80 1.9900634212544457
90 1.9866745407037671
100 1.983971518523551
125 1.9791241094237972
150 1.975905330896619
175 1.973612461954384
200 1.97189622363391
250 1.9694983934211534
300 1.967903011261086
350 1.9667650028636563
400 1.9659123432294678
450 1.965249664736467
500 1.9647198374673704
=#