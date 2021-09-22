using Distributions, Random
using LinearAlgebra
using DataFrames
using Plots
using RCall
using CSV
include("funcs.jl")

data = DataFrame(
        m20 = [3, 11, 10, 7, 4, 3, 2, 1, 0, 0, 0, 1, 1, 101],
        f20 = [0, 2, 4, 8, 9, 3, 0, 0, 0, 0, 0, 0, 0, 126],
        m32 = [7, 10, 11, 16, 3, 2, 1, 0, 0, 0, 0, 0, 0, 19],
        f32 = [1, 5, 11, 10, 5, 1, 0, 1, 0, 0, 0, 0, 0, 47],
        m50 = [5, 8, 11, 15, 4, 2, 1, 1, 0, 0, 0, 0, 0, 7],
        f50 = [0, 4, 6, 6, 3, 1, 1, 4, 0, 0, 0, 1, 1, 17],
        m80 = [4, 10, 8, 14, 8, 2, 1, 0, 0, 1, 0, 0, 0, 2],
        f80 = [2, 7, 15, 9, 3, 4, 1, 1, 0, 1, 0, 0, 0, 4]
)


β_p = Exponential(0.25)
α_p = Uniform(0,1)
γ_p = Gamma(1.5,1)
p_p = Gamma(1.5,1)

α_i = rand(α_p)
β_i = 0.00000001
γ_i = rand(γ_p)
p_i = rand(p_p)

# test of simple model
test = beetle_sim(144, 0.6, 5, 2, 2, 0.2)
println(test)
println(sum(test))

function test_prop()
    tot = zeros(14)
    for i = 1:10000
        tot += beetle_sim(144, 0.5, 5, 0.08, 0.31, 0.2)
    end
    return tot ./ 10000
end
# proportions from test of simple model
test_proportions = test_prop()

α_i = rand(α_p)
p_i = rand(p_p)
# test of extended model
ext_test = beetle_sim_ext(144, α_i, β_i, γ_i, p_i, 0.2)
println(ext_test)
println(sum(ext_test))


diggle = [3, 11, 10, 7, 4, 3, 2, 1, 0,0, 0, 1, 1, 101]

# abc simulation
sim = abc_sim(diggle, 0.2, 10^5, 144, α_p, β_p, γ_p, p_p)
