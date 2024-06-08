using Distributions, Random, DataFrames, CSV, StatsBase

# TODO replace with zscore from StatsBase
function standardize(x, μ, σ, n::Int64)::Float64
    (x - μ) / (σ / sqrt(n))
end

function standardize(data::Vector{Float64})::Vector{Float64}
    (data .- mean(data)) ./ std(data)
end

function t_score(data::Vector{Float64}; μ::Real)::Float64
    standardize(mean(data), μ, std(data), length(data))
end

function sampling_distribution(statistic::Function, d::Distribution, n::Int, r::Int, sample_statistics::Vector{Float64}; args...)::Vector{Float64}
    Random.seed!(0)
    # Preallocating sample vector for speed
    sample = zeros(n)

    # Sampling r times and calculating the statistic
    @inbounds for i in 1:r
        rand!(d, sample) # in-place to reduce memory allocation
        sample_statistics[i] = statistic(sample; args...)::Float64
    end

    # samples = zeros(n, r)
    # sample_statistics = mapcols(col -> statistic(col; args...), samples)

    sample_statistics
end

function analysis(statistic::Function, d::Distribution, n::Int, r::Int, μ::Real, σ::Real, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64}
    sample_statistics = zeros(r)
    statistics = sampling_distribution(statistic, d, n, r, sample_statistics; args...)
    skewness = StatsBase.skewness(statistics)
    kurtosis = StatsBase.kurtosis(statistics)

    # println("Original: ", minimum(statistics), " ", maximum(statistics))

    # Standardizing the values to look at tail probabilities
    # TODO fix standardization
    # z_scores = standardize.(statistics, μ, σ, 1)
    # println("Normalized:", minimum(z_scores), " ", maximum(z_scores))
    # println()

    # dt = fit(ZScoreTransform, z_scores, dims=1)
    # StatsBase.transform(dt, statistics)
    # println("Normalized Base: ", minimum(statistics), " ", maximum(statistics))
    # println()

    z_scores = standardize(statistics)
    # println("Normalized Base Other: ", minimum(statistics), " ", maximum(statistics))
    # println()

    upper = sum(z_scores .>= critical) / r
    lower = sum(z_scores .<= -critical) / r

    (upper, lower, skewness, kurtosis)
end

function analyze_distributions(statistic::Function, r::Int, sample_sizes::Vector{Int}, critical::Function, distributions, params=false)::DataFrame
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

    # Compile
    analyze_distributions(mean, 1, sample_sizes, zstar, distributions)
    analyze_distributions(t_score, 1, sample_sizes, tstar, distributions, true)

    # Warning: this code will take a very long time to run if used with a large r. We used r = 10_000_000
    @time means::DataFrame = analyze_distributions(mean, r, sample_sizes, zstar, distributions)
    CSV.write("means.csv", means)

    @time t::DataFrame = analyze_distributions(t_score, r, sample_sizes, tstar, distributions, true)
    CSV.write("t.csv", t)
end

function graphing(r)
    d = Exponential()
    μ = mean(d)
    σ = std(d)
    sample_statistics = zeros(r)
    n = 30

    # Creating sampling distributions
    exponential30 = sampling_distribution(statistic, d, n, r, sample_statistics)
    exponential30std = standardize.(exponential30, μ, σ, n)

    n = 150
    exponential150 = sampling_distribution(statistic, d, n, r, sample_statistics)
    exponential150std = standardize.(exponential150, μ, σ, n)

    graphing = DataFrame(
        "Exponential" => exponential30,
        "Exponential Z-Scores" => exponential30std,
        "Exponential 150" => exponential150,
        "Exponential 150 Z-Scores" => exponential150std
    )

    CSV.write("graphing1.csv", graphing) #todo change
end

main(10_000_000)
# graphing(10_000_000)