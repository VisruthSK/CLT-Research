using Distributions, Random, DataFrames, CSV, StatsBase

function t_statistic(data, μ, σ)
    (mean(data) - μ) / (σ / sqrt(length(data)))
end

function analysis(statistic::Function, d::Distribution, n::Int64, r::Int64, critical::Float64, μ::Float64, σ::Float64)
    Random.seed!(0)
    sample_statistics = zeros(r)
    sample = zeros(n)

    # TODO test again with this being threaded and locked
    @inbounds for i in 1:r
        rand!(d, sample)
        if statistic == mean
            sample_statistics[i] = statistic(sample)
        elseif statistic == t_statistic
            sample_statistics[i] = statistic(sample, μ, σ)
        end
    end

    m = mean(sample_statistics)
    s = std(sample_statistics)
    skewness = StatsBase.skewness(sample_statistics)
    kurtosis = StatsBase.kurtosis(sample_statistics)

    upper = sum(sample_statistics .>= m + critical * s) / r
    lower = sum(sample_statistics .<= m - critical * s) / r

    (upper, lower, m, s, skewness, kurtosis)
end

function analyze_distributions(statistic::Function, critical::Float64, r::Number)::DataFrame
    println("Analyzing distributions with $(r) repetitions")

    # sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]
    sample_sizes = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500]
    distributions = [
        LogNormal(0, 1.4865),
        Poisson(0.001),
        Gamma(0.02),
        LogNormal(0, 1),
        Exponential(),
        LogNormal(0, 0.5),
        LogNormal(0, 0.25),
        Beta(0.3, 0.2),
        Normal()
    ]
    results = DataFrame(
        Distribution=String[],
        Skewness=Float64[],
        Sample_Size=Int[],
        Upper_Tail=Float64[],
        Lower_Tail=Float64[],
        Tail_Sum=Float64[],
        Tail_Difference=Float64[],
        Sampling_Mean=Float64[],
        Sampling_SD=Float64[],
        Sampling_Skewness=Float64[],
        Sampling_Kurtosis=Float64[],
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
            upper, lower, m, s, sample_skew, sample_kurt = analysis(statistic, d, n, r, critical, μ, σ)
            Threads.lock(u) do
                push!(results, (string(d), skewness, n, upper, lower, upper + lower, upper - lower, m, s, sample_skew, sample_kurt, μ, σ))
            end
        end
    end

    sort!(results, [:Distribution, :Sample_Size])
end

# compile
analyze_distributions(mean, quantile(Normal(), 0.975), 1)
# analyze_distributions(t_statistic, quantile(Normal(), 0.975), 1)


@time CSV.write("test.csv", analyze_distributions(mean, quantile(Normal(), 0.975), 100000))

@time results = analyze_distributions(mean, quantile(Normal(), 0.975), 100000)
CSV.write("CLT_means.csv", results)

# @time results = analyze_distributions(t_statistic, quantile(TDist(100000 - 1), 0.975), 100000)
# CSV.write("CLT_t.csv", results)

# @code_warntype analyze_distributions(mean, quantile(Normal(), 0.975), 1000)

# using BenchmarkTools
# @profview analyze_distributions(mean, quantile(Normal(), 0.975), 10000)
# @btime analyze_distributions(mean, quantile(Normal(), 0.975), 1000)

# TODO Current implementation is not entirely reproducible, because of rand! and multithreading.