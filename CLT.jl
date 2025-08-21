using Distributions, Random, DataFrames, CSV, StatsBase, Statistics

function sampling_distribution(statistic::Function, d::Distribution, n::Int, r::Int; args...)::Vector{Float64}
    Random.seed!(0)
    # Preallocating vectors for speed
    sample_statistics = zeros(r)
    sample = zeros(n)

    # Sampling r times and calculating the statistic
    @inbounds for i in 1:r
        rand!(d, sample) # in-place to reduce memory allocation
        sample_statistics[i] = statistic(sample; args...)::Float64
    end

    sample_statistics
end

function analysis(statistic::Function, d::Distribution, n::Int, r::Int, μ::Real, σ::Real, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    statistics = sampling_distribution(statistic, d, n, r; args...)
    skewness = StatsBase.skewness(statistics)
    kurtosis = StatsBase.kurtosis(statistics)

    # NOT STANDARDIZED
    md = median(statistics)
    IQR = quantile(statistics, 0.75) - quantile(statistics, 0.25)

    # Standardizing the values to look at tail probabilities
    z_scores = zeros(r)
    zscore!(z_scores, statistics, μ, σ / sqrt(n))

    # Calculating tail probabilities
    upper = sum(z_scores .>= critical) / r
    lower = sum(z_scores .<= -critical) / r

    (md, IQR, upper, lower, skewness, kurtosis)
end

function analyze_distributions(statistic::Function, r::Int, sample_sizes::Vector{Int}, critical::Function, distributions, params=false)::DataFrame
    println("Analyzing sampling distributions of $(statistic)s with $(r) repetitions")
    # Setting up the results we're interested in
    results = DataFrame(
        "Distribution" => String[],
        "Skewness" => Float64[],
        "Sample Size" => Int64[],
        "Median" => Float64[],
        "IQR" => Float64[],
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
                median, IQR, upper, lower, sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)), μ=μ)
            else
                median, IQR, upper, lower, sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)))
            end
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, median, IQR, upper, lower, sample_skewness, sample_kurtosis, μ, σ))
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
    # Compile
    analyze_distributions(mean, 1, sample_sizes, zstar, distributions)

    means::DataFrame = analyze_distributions(mean, r, sample_sizes, quantile(Normal(), 0.975), distributions)
    CSV.write("means.csv", means, compress=true)

    nothing
end

function graphing(r)
    d = Exponential()
    μ = mean(d)
    σ = std(d)

    n = 30
    exponential30 = sampling_distribution(mean, d, n, r)
    exponential30std = zeros(r)
    zscore!(exponential30std, exponential30, μ, σ / sqrt(n))

    n = 150
    exponential150 = sampling_distribution(mean, d, n, r)
    exponential150std = zeros(r)
    zscore!(exponential150std, exponential150, μ, σ / sqrt(n))

    d = Normal()
    μ = mean(d)
    σ = std(d)

    n = 5
    normal5 = sampling_distribution(mean, d, n, r)
    normal5std = zeros(r)
    zscore!(normal5std, normal5, μ, σ / sqrt(n))

    n = 10
    normal10 = sampling_distribution(mean, d, n, r)
    normal10std = zeros(r)
    zscore!(normal10std, normal10, μ, σ / sqrt(n))

    n = 30
    normal30 = sampling_distribution(mean, d, n, r)
    normal30std = zeros(r)
    zscore!(normal30std, normal30, μ, σ / sqrt(n))

    graphing = DataFrame(
        "Exponential 30" => exponential30,
        "Exponential 30 Z-Scores" => exponential30std,
        "Exponential 150" => exponential150,
        "Exponential 150 Z-Scores" => exponential150std,
        "Normal 5" => normal5,
        "Normal 5 Z-Scores" => normal5std,
        "Normal 10" => normal10,
        "Normal 10 Z-Scores" => normal10std,
        "Normal 30" => normal30,
        "Normal 30 Z-Scores" => normal30std,
    )

    CSV.write("graphing.csv", graphing, compress=true)

    nothing
end

function gamma_graphing(r)
    d = Gamma(16, 1)
    μ = mean(d)
    σ = std(d)

    n = 10
    gamma10 = sampling_distribution(mean, d, n, r)
    gamma10std = zeros(r)
    zscore!(gamma10std, gamma10, μ, σ / sqrt(n))

    graphing = DataFrame(
        "Gamma 10" => gamma10,
        "Gamma 10 Z-Scores" => gamma10std,
    )

    CSV.write("gamma_graphing.csv", graphing, compress=true)

    nothing
end

main(1_000_000)
graphing(1_000_000)
gamma_graphing(1_000_000)
