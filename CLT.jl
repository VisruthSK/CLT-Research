using Distributions, Random, DataFrames, CSV

function t_statistic(data, μ, σ)
    n = length(data)
    (mean(data) - μ) / (σ / sqrt(n))
end

function analysis(statistic, distro, n, r, params...)
    Random.seed!(0)

    d = distro(params...)
    μ = mean(d)
    σ = std(d)

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

    # TODO not extensible
    critical = statistic == mean ? quantile(Normal(), 0.975) : quantile(TDist(n - 1), 0.975)

    upper = sum(sample_statistics .>= m + critical * s) / r
    lower = sum(sample_statistics .<= m - critical * s) / r

    (upper, lower, upper + lower, upper - lower)
end

function analyze_distributions(statistic, r)
    println("Analyzing distributions with $(r) repetitions")

    sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]
    distributions = [
        (LogNormal, [0, 1.4865], 31.65266),
        (Poisson, [0.001], 31.6228),
        (Gamma, [0.02], 20),
        (LogNormal, [0, 1], 6.1849),
        (Exponential, [], 2),
        (LogNormal, [0, 0.5], 1.7502),
        (LogNormal, [0, 0.25], 0.7883),
        (Beta, [0.3, 0.2], -0.4),
        (Normal, [], 0)
    ]
    results = DataFrame(Distribution=String[], Params=Vector[], Skewness=Float64[], Sample_Size=Int[], Upper_Tail=Float64[], Lower_Tail=Float64[], Tail_Sum=Float64[], Tail_Difference=Float64[])

    for (distro, params, skewness) in distributions
        println("$(distro) with parameters $(params)")
        for n in sample_sizes
            upper, lower, tail_sum, tail_diff = analysis(statistic, distro, n, r, params...)
            push!(results, (string(distro), params, skewness, n, upper, lower, tail_sum, tail_diff))
        end
    end

    results
end

analyze_distributions(mean, 1) # compile

results = analyze_distributions(mean, 100000)
CSV.write("CLT_means.csv", results)

results = analyze_distributions(t_statistic, 100000)
CSV.write("CLT_t.csv", results)