using Distributions, Random, DataFrames, CSV, StatsBase

function t_statistic(data::Vector, μ::Float64)
    (mean(data) - μ) / (std(data) / sqrt(length(data)))
end

function analysis(statistic::Function, d::Distribution, n::Int64, r::Int64, μ::Number)
    Random.seed!(0)
    sample_statistics = zeros(r)
    sample = zeros(n)

    @inbounds for i in 1:r
        rand!(d, sample)
        if statistic == mean
            sample_statistics[i] = statistic(sample)
        elseif statistic == t_statistic
            sample_statistics[i] = statistic(sample, μ)
        end
    end

    skewness = StatsBase.skewness(sample_statistics)

    μ::Float64 = mean(d)
    σ::Float64 = std(d)

    z_scores = (sample_statistics .- μ) ./ (σ / sqrt(n))
    upper = sum(z_scores .>= 1.96) / r
    lower = sum(z_scores .<= -1.96) / r

    (upper, lower, skewness)
end

function analyze_distributions(statistic::Function, r::Number)::DataFrame
    println("Analyzing distributions with $(r) repetitions")

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

    @inbounds for d::Distribution in distributions
        println("$(string(d))")

        μ::Float64 = mean(d)
        σ::Float64 = std(d)
        skewness::Float64 = StatsBase.skewness(d)

        u = Threads.SpinLock()
        @inbounds Threads.@threads for n in sample_sizes
            upper, lower, sample_skew = analysis(statistic, d, n, r, μ)
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, upper, lower, sample_skew, μ, σ))
            end
        end

    end

    sort!(results, [:Distribution, :Sample_Size])
end

# compile
analyze_distributions(mean, 1)

results = analyze_distributions(mean, 1000000)
CSV.write("means.csv", results)