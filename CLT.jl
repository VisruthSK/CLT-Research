using Distributions, Random, DataFrames, CSV, StatsBase, Statistics

function t_statistic(data::Vector{Float64}; μ::Real)::Float64
    (mean(data) - μ) / (std(data) / sqrt(length(data)))
end

function sampling_distribution(statistic::Function, d::Distribution, n::Int, r::Int; args...)::Tuple{Vector{Float64},Float64}
    rng = MersenneTwister(0)
    # Preallocating vectors for speed
    sample_statistics = zeros(r)
    sample = zeros(n)
    sample_skewness_total = 0.0

    # Sampling r times and calculating the statistic
    @inbounds for i in 1:r
        rand!(rng, d, sample) # in-place to reduce memory allocation
        sample_statistics[i] = statistic(sample; args...)::Float64
        sample_skewness_total += StatsBase.skewness(sample)
    end

    (sample_statistics, sample_skewness_total / r)
end

function analysis(statistic::Function, d::Distribution, n::Int, r::Int, μ::Real, σ::Real, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64}
    statistics, average_sample_skewness = sampling_distribution(statistic, d, n, r; args...)
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

    (md, IQR, upper, lower, skewness, average_sample_skewness, kurtosis)
end

function analysis_studentized(statistic::Function, d::Distribution, n::Int, r::Int, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64}
    statistics, average_sample_skewness = sampling_distribution(statistic, d, n, r; args...)
    skewness = StatsBase.skewness(statistics)
    kurtosis = StatsBase.kurtosis(statistics)

    md = median(statistics)
    IQR = quantile(statistics, 0.75) - quantile(statistics, 0.25)

    upper = sum(statistics .>= critical) / r
    lower = sum(statistics .<= -critical) / r

    (md, IQR, upper, lower, skewness, average_sample_skewness, kurtosis)
end

function sampling_distribution_two_groups(statistic::Function, d::Distribution, n1::Int, n2::Int, r::Int; args...)::Tuple{Vector{Float64},Float64}
    rng = MersenneTwister(0)
    difference_statistics = zeros(r)
    sample_1 = zeros(n1)
    sample_2 = zeros(n2)
    sample_skewness_total = 0.0

    @inbounds for i in 1:r
        rand!(rng, d, sample_1)
        rand!(rng, d, sample_2)
        g1 = statistic(sample_1; args...)::Float64
        g2 = statistic(sample_2; args...)::Float64
        difference_statistics[i] = g1 - g2
        sample_skewness_total += (StatsBase.skewness(sample_1) + StatsBase.skewness(sample_2)) / 2
    end

    (difference_statistics, sample_skewness_total / r)
end

function analysis_difference_in_means(statistic::Function, d::Distribution, n1::Int, n2::Int, r::Int, σ::Real, critical::Float64; args...)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64}
    difference_statistics, average_sample_skewness = sampling_distribution_two_groups(statistic, d, n1, n2, r; args...)
    skewness = StatsBase.skewness(difference_statistics)
    kurtosis = StatsBase.kurtosis(difference_statistics)

    md = median(difference_statistics)
    IQR = quantile(difference_statistics, 0.75) - quantile(difference_statistics, 0.25)

    z_scores = zeros(r)
    sd_difference = σ * sqrt(1 / n1 + 1 / n2)
    zscore!(z_scores, difference_statistics, 0.0, sd_difference)

    upper = sum(z_scores .>= critical) / r
    lower = sum(z_scores .<= -critical) / r

    (md, IQR, upper, lower, skewness, average_sample_skewness, kurtosis)
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
        "Average Sample Skewness" => Float64[],
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
                median, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)), μ=μ)
            else
                median, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ, abs(critical(n)))
            end
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, median, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis, μ, σ))
            end
        end

    end

    sort!(results, [:Distribution, :"Sample Size"])
end

function analyze_t_statistics(r::Int, sample_sizes::Vector{Int}, critical::Function, distributions)::DataFrame
    println("Analyzing sampling distributions of t statistics with $(r) repetitions")
    results = DataFrame(
        "Distribution" => String[],
        "Skewness" => Float64[],
        "Sample Size" => Int64[],
        "Median" => Float64[],
        "IQR" => Float64[],
        "Upper Tail" => Float64[],
        "Lower Tail" => Float64[],
        "Sampling Skewness" => Float64[],
        "Average Sample Skewness" => Float64[],
        "Sampling Kurtosis" => Float64[],
        "Population Mean" => Float64[],
        "Population SD" => Float64[]
    )

    @inbounds for d::Distribution in distributions
        println(string(d))

        μ = mean(d)
        σ = std(d)
        skewness = StatsBase.skewness(d)

        u = Threads.SpinLock()
        @inbounds Threads.@threads for n in sample_sizes
            median, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis = analysis_studentized(
                t_statistic, d, n, r, abs(critical(n)), μ=μ
            )
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, median, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis, μ, σ))
            end
        end
    end

    sort!(results, [:Distribution, :"Sample Size"])
end

function analyze_difference_in_means(statistic::Function, r::Int, sample_sizes::Vector{Int}, ratios::Vector{Tuple{Int,Int}}, critical::Function, distributions)::DataFrame
    println("Analyzing sampling distributions of difference in $(statistic)s with $(r) repetitions")
    summary = DataFrame(
        "Distribution" => String[],
        "Skewness" => Float64[],
        "Sample Size 1" => Int64[],
        "Sample Size 2" => Int64[],
        "Median" => Float64[],
        "IQR" => Float64[],
        "Upper Tail" => Float64[],
        "Lower Tail" => Float64[],
        "Sampling Skewness" => Float64[],
        "Average Sample Skewness" => Float64[],
        "Sampling Kurtosis" => Float64[],
        "Population Mean" => Float64[],
        "Population SD" => Float64[]
    )

    for d::Distribution in distributions
        println(string(d))
        μ = mean(d)
        σ = std(d)
        skewness = StatsBase.skewness(d)

        for n in sample_sizes
            for ratio in ratios
                n1 = n
                n2 = cld(ratio[1] * n, ratio[2])

                md, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis = analysis_difference_in_means(
                    statistic, d, n1, n2, r, σ, abs(critical(n1, n2))
                )
                push!(
                    summary,
                    (string(d), skewness, n1, n2, md, IQR, upper, lower, sampling_skewness, average_sample_skewness, sample_kurtosis, μ, σ)
                )
            end
        end
    end

    sort!(summary, [:Distribution, :"Sample Size 1", :"Sample Size 2"])
    summary
end

function main(r; mode::Symbol=:means, sample_sizes::Vector{Int}=[5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500])
    distributions = [
        Gamma(16),
        LogNormal(0, 0.25),
        Gamma(4),
        Gamma(2),
        LogNormal(0, 0.5),
        # Gamma(1),
        Exponential(),
        Gamma(0.64),
        LogNormal(0, 0.75)
    ]
    ratios = [(1, 1), (1, 2), (1, 3)]
    zstar = quantile(Normal(), 0.975)
    tstar = n -> quantile(TDist(n - 1), 0.975)
    critical_mean = n -> zstar
    critical_diff = (n1, n2) -> zstar
    valid_modes = (:means, :difference_in_means, :both, :t, :all)

    mode in valid_modes || throw(ArgumentError("Unsupported mode $(mode). Valid modes are $(valid_modes)."))

    if mode in (:means, :both, :all)
        means::DataFrame = analyze_distributions(mean, r, sample_sizes, critical_mean, distributions)
        CSV.write("means.csv", means, compress=true)
    end

    if mode in (:difference_in_means, :both, :all)
        difference_summary = analyze_difference_in_means(mean, r, sample_sizes, ratios, critical_diff, distributions)
        CSV.write("difference_means.csv", difference_summary, compress=true)
    end

    if mode in (:t, :all)
        t_summary::DataFrame = analyze_t_statistics(r, sample_sizes, tstar, distributions)
        CSV.write("t_statistics.csv.gz", t_summary, compress=true)
    end

    nothing
end

function graphing(r)
    d = Exponential()
    μ = mean(d)
    σ = std(d)

    n = 30
    exponential30, _ = sampling_distribution(mean, d, n, r)
    exponential30std = zeros(r)
    zscore!(exponential30std, exponential30, μ, σ / sqrt(n))

    n = 150
    exponential150, _ = sampling_distribution(mean, d, n, r)
    exponential150std = zeros(r)
    zscore!(exponential150std, exponential150, μ, σ / sqrt(n))

    d = Normal()
    μ = mean(d)
    σ = std(d)

    n = 5
    normal5, _ = sampling_distribution(mean, d, n, r)
    normal5std = zeros(r)
    zscore!(normal5std, normal5, μ, σ / sqrt(n))

    n = 10
    normal10, _ = sampling_distribution(mean, d, n, r)
    normal10std = zeros(r)
    zscore!(normal10std, normal10, μ, σ / sqrt(n))

    n = 30
    normal30, _ = sampling_distribution(mean, d, n, r)
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
    gamma10, _ = sampling_distribution(mean, d, n, r)
    gamma10std = zeros(r)
    zscore!(gamma10std, gamma10, μ, σ / sqrt(n))

    graphing = DataFrame(
        "Gamma 10" => gamma10,
        "Gamma 10 Z-Scores" => gamma10std,
    )

    CSV.write("gamma_graphing.csv", graphing, compress=true)

    nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # main(1_000_000; mode=:means)
    main(1_000_000; mode=:t, sample_sizes=[5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000])
    # main(1_000_000; mode=:difference_in_means)
    # graphing(1_000_000)
    # gamma_graphing(1_000_000)
end
