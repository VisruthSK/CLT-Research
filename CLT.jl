using Distributions, Random, DataFrames, CSV, StatsBase

function t_statistic(data, μ, σ)
    n = length(data)
    (mean(data) - μ) / (σ / sqrt(n))
end

function analysis(statistic, d, n, r, critical, μ, σ)
    Random.seed!(0)

    sample_statistics = zeros(r)
    Threads.@threads for i in 1:r
        sample_statistics[i] = try
            statistic(rand(d, n), μ, σ)
        catch
            statistic(rand(d, n))
        end
    end

    m = mean(sample_statistics)
    s = std(sample_statistics)

    c = critical
    upper = sum(sample_statistics .>= m + c * s) / r
    lower = sum(sample_statistics .<= m - c * s) / r

    (upper, lower, upper + lower, upper - lower, m, s, μ, σ)
end

function analyze_distributions(statistic, critical, r)
    println("Analyzing distributions with $(r) repetitions")

    sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]
    # distributions = [
    #     (() -> LogNormal(0, 1.4865), 31.65266),
    #     (() -> Poisson(0.001), 31.6228),
    #     (() -> Gamma(0.02), 20),
    #     (() -> LogNormal(0, 1), 6.1849),
    #     (() -> Exponential(), 2),
    #     (() -> LogNormal(0, 0.5), 1.7502),
    #     (() -> LogNormal(0, 0.25), 0.7883),
    #     (() -> Beta(0.3, 0.2), -0.4),
    #     (() -> Normal(), 0)
    # ]

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
    results = DataFrame(Distribution=String[], Skewness=Float64[], Sample_Size=Int[], Upper_Tail=Float64[], Lower_Tail=Float64[], Tail_Sum=Float64[], Tail_Difference=Float64[], Sampling_Mean=Float64[], sampling_SD=Float64[], Population_Mean=Float64[], Population_SD=Float64[])

    for d in distributions
        println("$(string(d))")
        # d = distro_gen()
        μ = mean(d)
        σ = std(d)

        skewness = StatsBase.skewness(d)
        for n in sample_sizes
            upper, lower, tail_sum, tail_diff, m, s, μ, σ = analysis(statistic, d, n, r, critical, μ, σ)
            push!(results, (string(d), skewness, n, upper, lower, tail_sum, tail_diff, m, s, μ, σ))
        end
    end

    results
end

analyze_distributions(mean, quantile(Normal(), 0.975), 1) # compile
analyze_distributions(t_statistic, quantile(Normal(), 0.975), 1) # compile

@time results = analyze_distributions(mean, quantile(Normal(), 0.975), 100000)
CSV.write("CLT_means.csv", results)

@time results = analyze_distributions(t_statistic, quantile(TDist(100000 - 1), 0.975), 100000)
CSV.write("CLT_t.csv", results)