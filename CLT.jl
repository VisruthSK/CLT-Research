using Distributions, Random, DataFrames, CSV, StatsBase

function t_statistic(data::Vector, μ::Float64)
    (mean(data) - μ) / (std(data) / sqrt(length(data)))
end

function analysis(statistic::Function, d::Distribution, n::Int64, r::Int64, μ::Float64, σ::Float64)::{Float64,Float64,Float64}
    # Setting some basic variables
    Random.seed!(0)
    sample_statistics = zeros(r)
    sample = zeros(n)

    # Create r samples of size n and store their statistics
    @inbounds for i in 1:r
        rand!(d, sample)
        sample_statistics[i] = statistic(sample)
    end

    # Getting the skewness of the sampling distribution
    skewness = StatsBase.skewness(sample_statistics)

    # Standardizing and calculating the proportion of sample statistics that fall in the tails
    z_scores = (sample_statistics .- μ) ./ (σ ./ sqrt(n))
    upper = sum(z_scores .>= 1.96) / r
    lower = sum(z_scores .<= -1.96) / r

    (upper, lower, skewness)
end

function analyze_distributions(statistic::Function, r::Int64)::DataFrame
    println("Analyzing distributions with $(r) repetitions")

    # Some setup
    sample_sizes = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500]
    distributions = [
        Gamma(16),
        LogNormal(0, 0.25),
        Gamma(4),
        Gamma(2),
        LogNormal(0, 0.5),
        Gamma(1),
        Exponential(),
        LogNormal(0, 0.75)
    ]
    results = DataFrame(
        Distribution=String[],
        Skewness=Float64[],
        Sample_Size=Int64[],
        Upper_Tail=Float64[],
        Lower_Tail=Float64[],
        Sampling_Skewness=Float64[],
        Population_Mean=Float64[],
        Population_SD=Float64[]
    )

    # Analyzing each distribution
    @inbounds for d::Distribution in distributions
        println("$(string(d))")

        # Getting population parameters
        μ = mean(d)
        σ = std(d)
        skewness = StatsBase.skewness(d)

        u = Threads.SpinLock() # lock to avoid data races
        @inbounds Threads.@threads for n in sample_sizes
            upper, lower, sample_skew = analysis(statistic, d, n, r, μ, σ)
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, upper, lower, sample_skew, μ, σ))
            end
        end

    end

    sort!(results, [:Distribution, :Sample_Size])
end

function main()
    # Compile
    analyze_distributions(mean, 1)

    # ~50 seconds on a i7-12700h
    results = analyze_distributions(mean, 1_000_000)
    CSV.write("means.csv", results)
end

main()