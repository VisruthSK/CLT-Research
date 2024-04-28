using Distributions, Random, DataFrames, TidierPlots

# TODO Change to 100k
const r::Number = 1000
const sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]

function analysis(statistic, distro, n, params...)
    Random.seed!(0)
    d = distro(params...)

    sample_statistics = zeros(r)

    # sample_statistics = pmap(x -> statistic(rand(d, n)), 1:r)

    # sample_statistics = [statistic(rand(d, n)) for i in 1:r]

    Threads.@threads for i in 1:r
        sample_statistics[i] = statistic(rand(d, n))
    end

    s = std(sample_statistics)
    m = mean(sample_statistics)

    upper = sum(sample_statistics .>= m + 1.96 * s) / r
    lower = sum(sample_statistics .<= m - 1.96 * s) / r

    (upper, lower, upper + lower, upper - lower)
end

const distributions = [
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

# time this block
@time begin
    for (distro, params, skewness) in distributions
        println("$(distro) with parameters $(params)")
        for n in sample_sizes
            upper, lower, tail_sum, tail_diff = analysis(mean, distro, n, params...)
            push!(results, (string(distro), params, skewness, n, upper, lower, tail_sum, tail_diff))
        end
    end
end

display(results)

filtered = filter(row -> row[:Params] in [[]], results)

ggplot(filtered, aes(x=:Sample_Size, y=:Tail_Sum, color=:Distribution)) +
geom_line() +
geom_point() +
theme_minimal() +
labs(title="Central Limit Theorem", x="Sample Size", y="Tail Sum")

# TODO Format results in table
# TODO Plot results, lin reg?