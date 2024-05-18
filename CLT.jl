using Distributions, Random, DataFrames, CSV, StatsBase


# TODO replace with zscore from StatsBase
function standardize(x, μ, σ, n::Int64)::Float64
    (x - μ) / (σ / sqrt(n))
end

function t_score(data::Vector{Float64}, μ::Float64)::Float64
    standardize(mean(data), μ, std(data), length(data))
end

function sampling_distribution(d::Distribution, n::Int, r::Int; statistic::Function, args...)::Vector{Float64}
    Random.seed!(0)
    # Preallocating vectors for speed
    sample_statistics = zeros(r)
    sample = zeros(n)

    # Sampling r times and calculating the statistic
    @inbounds for i in 1:r
        rand!(d, sample) # in-place to reduce memory allocation

        # sample_statistics[i] = statistic(sample, args...)::Float64

        # TODO Hacky, fix
        try
            # Threads.threadid() == 1 ? println(args) : nothing
            sample_statistics[i] = statistic(sample; args)::Float64
        catch e
            if isa(e, MethodError)
                sample_statistics[i] = statistic(sample)::Float64
            else
                rethrow(e)
            end
        end
    end

    sample_statistics
end

function analysis(statistic, d::Distribution, n::Int, r::Int, μ::Real, σ::Real)::Tuple{Float64,Float64,Float64,Float64}
    sample_statistics = sampling_distribution(d, n, r, statistic=statistic, μ=μ)
    skewness = StatsBase.skewness(sample_statistics)
    kurtosis = StatsBase.kurtosis(sample_statistics)

    # Standardizing the values to look at tail probabilities
    z_scores = standardize.(sample_statistics, μ, σ, n)
    # TODO change 1.96 to variable
    upper = mean(z_scores .>= 1.96)
    lower = mean(z_scores .<= 1.96)
    # lower = sum(z_scores .<= -1.96) / r

    (upper, lower, skewness, kurtosis)
end

function analyze_distributions(statistic, r::Int64, sample_sizes::Vector{Int64}, distributions)::DataFrame
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
            upper, lower, sample_skewness, sample_kurtosis = analysis(statistic, d, n, r, μ, σ)
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, upper, lower, sample_skewness, sample_kurtosis, μ, σ))
            end
        end

    end

    sort!(results, [:Distribution, :"Sample Size"])
end

function main(r)
    sample_sizes = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500]
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
    analyze_distributions(mean, 1, sample_sizes, distributions)
    analyze_distributions(t_score, 1, sample_sizes, distributions)

    # Warning: this code will take a very long time to run if used with a large r. We used r = 10_000_000
    means::DataFrame = analyze_distributions(mean, r, sample_sizes, distributions)
    CSV.write("means.csv", means)

    t::DataFrame = analyze_distributions(t_score, r, sample_sizes, distributions)
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

main(10_000_000)
graphing(10_000_000)
