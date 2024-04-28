using Distributions, Random, DataFrames, CSV

function analysis(statistic, distro, n, r, params...)
    Random.seed!(0)
    d = distro(params...)

    sample_statistics = zeros(r)
    Threads.@threads for i in 1:r
        sample_statistics[i] = statistic(rand(d, n))
    end

    m = mean(sample_statistics)
    s = std(sample_statistics)

    upper = sum(sample_statistics .>= m + 1.96 * s) / r
    lower = sum(sample_statistics .<= m - 1.96 * s) / r

    (upper, lower, upper + lower, upper - lower)
end

function analyze_distributions(r)
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
            upper, lower, tail_sum, tail_diff = analysis(mean, distro, n, r, params...)
            push!(results, (string(distro), params, skewness, n, upper, lower, tail_sum, tail_diff))
        end
    end

    results
end

analyze_distributions(1) # compile
results = analyze_distributions(10000)
display(results)

CSV.write("CLT.csv", results)