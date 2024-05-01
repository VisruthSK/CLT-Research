using Distributions, Random, DataFrames, CSV, StatsBase

function t_statistic(data, μ, σ)
    (mean(data) - μ) / (σ / sqrt(length(data)))
end

function analysis(statistic, d, n, r, critical, μ, σ)
    Random.seed!(0)

    sample_statistics = zeros(r)

    @inbounds Threads.@threads for i in 1:r
        rand!(d, sample_statistics)
        # sample = rand(d, n)
        # if statistic == mean
        #     sample_statistics[i] = statistic(sample)
        # elseif statistic == t_statistic
        #     sample_statistics[i] = statistic(sample, μ, σ)
        # end
    end

    if statistic == mean
        sample_statistics = statistic.(sample_statistics)
    elseif statistic == t_statistic
        sample_statistics = statistic.(sample_statistics, μ, σ)
    end

    m = mean(sample_statistics)
    s = std(sample_statistics)

    upper = sum(sample_statistics .>= m + critical * s) / r
    lower = sum(sample_statistics .<= m - critical * s) / r

    (upper, lower, m, s)
end

function analyze_distributions(statistic, critical, r)
    println("Analyzing distributions with $(r) repetitions")

    sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]

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

        μ = mean(d)
        σ = std(d)
        skewness = StatsBase.skewness(d)

        for n in sample_sizes
            upper, lower, m, s = analysis(statistic, d, n, r, critical, μ, σ)
            push!(results, (string(d), skewness, n, upper, lower, upper + lower, upper - lower, m, s, μ, σ))
        end
    end

    results
end

# compile
analyze_distributions(mean, quantile(Normal(), 0.975), 1)
analyze_distributions(t_statistic, quantile(Normal(), 0.975), 1)

@time results = analyze_distributions(mean, quantile(Normal(), 0.975), 100000)
CSV.write("CLT_means.csv", results)

@time results = analyze_distributions(t_statistic, quantile(TDist(100000 - 1), 0.975), 100000)
CSV.write("CLT_t.csv", results)



using BenchmarkTools
@profview analyze_distributions(mean, quantile(Normal(), 0.975), 1000)
@btime analyze_distributions(mean, quantile(Normal(), 0.975), 1000)