using Distributions, Random, DataFrames

# TODO Change to 100k
r::Number = 100
sample_sizes = [5, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000, 5000, 10000]

function analysis(statistic, distro, n, params...)
    Random.seed!(0)
    d = distro(params...)

    sample_statistics = zeros(r)
    # TODO see if threading makes sense here
    Threads.@threads for i in 1:r
        sample_statistics[i] = statistic(rand(d, n))
    end

    s = std(sample_statistics)
    m = mean(sample_statistics)

    upper = sum(sample_statistics .>= m + 1.96 * s) / r
    lower = sum(sample_statistics .<= m - 1.96 * s) / r

    upper - lower
end

distributions = [
    (LogNormal, [0, 1.4865]),
    (Poisson, [0.001]),
    (Gamma, [0.02]),
    (LogNormal, [0, 1]),
    (Exponential, []),
    (LogNormal, [0, 0.5]),
    (LogNormal, [0, 0.25]),
    (Beta, [0.3, 0.2]),
    (Normal, [])
]

results = DataFrame(Distribution=String[], Sample_Size=Int[], Tail_Difference=Float64[])

using BenchmarkTools
@btime begin
    for (distro, params) in distributions
        println("$(distro) with parameters $(params)")
        for n in sample_sizes
            tail_diff = analysis(mean, distro, n, params...)
            push!(results, ["$(distro) $(params)", n, tail_diff])
        end
    end
end


df_wide = unstack(results, :Sample_Size, :Distribution, :Tail_Difference)

results = combine(groupby(results, :Sample_Size), :Tail_Difference => (x -> mean(x)) => :Mean_Tail_Difference)

# Optionally, rename the columns to include the distribution names
rename!(results, Symbol("$(unique(results.Distribution))") .=> :Mean_Tail_Difference)

display(results)

# for (distro, params) in distributions
#     for n in sample_sizes
#         println("n = $n: $(analysis(mean, distro, n, params...))")
#     end
# end

# df <- data.frame(
#   sapply(sample_sizes, analysis, statistic = mean, distro = rlnorm, meanlog = 0, sdlog = 1.4865),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rpois, lambda = 0.001),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rgamma, shape = 0.02),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rlnorm, meanlog = 0, sdlog = 1),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rexp),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rlnorm, meanlog = 0, sdlog = 0.5),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rlnorm, meanlog = 0, sdlog = 0.25),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rbeta, shape1 = 0.3, shape2 = 0.2),
#   sapply(sample_sizes, analysis, statistic = mean, distro = rnorm)
# )