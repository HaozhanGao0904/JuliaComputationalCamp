############ StaticArrays ############

# arrays in Julia and most languages have a fixed number of dimension (type), but not a fixed size
M = rand(2,2)
M2 = rand(3,3)

typeof(M) == typeof(M2) # this will slow down Julia, because it doesn't know their sizes

# example
function add_up(M,R)
    total = 0.0
    
    for i = 1:1000000
        total += M[1,:]'*R
        total += M[2,:]'*R
    end
    return total
end

function add_up2(M,R)
    total = 0.0

    for i = 1:1000000
        total += M[1,1]*R[1] + M[1,2]*R[2]
        total += M[2,1]*R[1] + M[1,2]*R[2]
    end

    return total

end

using BenchmarkTools

M = rand(2,2)
R = rand(2)

@benchmark add_up(M,R)
@benchmark add_up2(M,R)

# we can fix this problem by using package StaticArrays
using StaticArrays
M_static = SMatrix{2,2}(M)
R_static = SVector{2}(R)

@benchmark add_up(M_static,R_static)

# StaticArrays is only good for "small" vectors and matrices. Rule of thumb is less than 100 numbers.

# Packages
using Random, Plots, Distributions, Statistics, Parameters

############ birthday problem: find the prob. that 2 ppl in a room with N people have same birthday ############

function birthday(N, sims)
    results = zeros(sims)

    for i = 1:sims
        # random birthdays
        days = rand(1:365, N)
        results[i] = length(unique(days))
    end

    results
    
end
# suppose there are 20 people in a room
res_20 = birthday(20, 10^6)
histogram(res_20)
# chance that 2 people have the same birthday
1.0 - mean(res_20.==20)

# suppose there are 70 people in a room
res_70 = birthday(70,10^6)
1- mean(res_70.==70)

############ average distance between 2 pts in a cube ############
 
function point_distance(sims)
    results = zeros(sims)

    for i = 1:sims
        p1 = rand(3)
        p2 = rand(3)
        
        results[i] = sqrt(sum((p1.-p2).^2))
    end
    return results  
end

# 1 million draws
res = point_distance(10^6)
mean(res)

# 100 draws
res_100 = point_distance(100)
mean(res_100)

############ expected value college given wage offer shock ############
@with_kw struct Primitives
    # w = exp(β_0 + β_1*s + ϵ)
    β_0::Float64 = 2.7
    β_1::Float64 = 0.47

    σ::Float64 = 0.597
    α::Float64 = 1.0
    B::Float64 = 5.0
    d::Float64 = 0.25
end

@with_kw struct Results
    emax::Vector{Float64}
    lfp::Vector{Float64} # labor force participation
    ewage::Vector{Float64} # expected wage
    ewage_obs::Vector{Float64}
end

# initialize model primitives and execute solutions
function Solve_Model(sims)
    prim = Primitives()
    res = Results(zeros(2), zeros(2), zeros(2), zeros(2))
    Compute_emax(prim, res, sims)
    return prim, res    
end

function Compute_emax(prim, res, sims)
    @unpack β_0, β_1, σ, α, B, d = prim
    dist = Normal(0.0, σ)

    val_notwork = α + log(B)
    
    # utility, lfp, wages based on schooling decision
    utils = zeros(2, sims)
    lfps = zeros(2, sims)
    wages = zeros(2, sims)

    # loop over possible schooling decisions
    for s = 0:1
        for i = 1:sims
            ϵ = rand(dist)
            wage = exp(β_0 + β_1 * s + ϵ)
            util = max(log(wage), val_notwork)
            utils[s+1, i] = util # update utility
            lfps[s+1,i] = (log(wage) > val_notwork) # updata labor force utility decision
            wages[s+1,i] = wage # update wage
        end
    end
    # expected continuation value
    res.emax[1] = mean(utils[1,:])
    res.emax[2] = mean(utils[2,:])
    # probability of working
    res.lfp[1] = mean(lfps[1,:])
    res.lfp[2] = mean(lfps[2,:])
    # average wage
    res.ewage[1] = mean(wages[1,:])
    res.ewage[2] = mean(wages[2,:])

    # wages of the people that we would see in data
    res.ewage_obs[1] = mean(wages[1,:].*lfps[1,:])/res.lfp[1]
    res.ewage_obs[2] = mean(wages[2,:].*lfps[2,:])/res.lfp[2]

end

prim, res = Solve_Model(100)

res.emax
res.lfp
res.ewage
res.ewage_obs # this is biased upwards because they chose to work

res.ewage[2] - res.ewage[1]
res.ewage_obs[2] - res.ewage_obs[1]
# college wage premium is biased downward because people who don't go to college are more likely to choose not to work

############# Quadrature #############

using FastGaussQuadrature, LinearAlgebra

# Gauss-Legendre Quadrature for approximating integrals between -1 and 1
f(x) = 4*x^3 - 9*x^2 - 8*x
F(x) = x^4 - 3*x^3 - 4*x^2

# calculate analytically
F(1) - F(-1)

# nodes and weights: Number of nodes and weights based on the degree of your polynomial
x, w = gausslegendre(3)
x
w
I = sum(w.*f.(x))


# Gauss-Hermite Quadrature for Normal Distributions
# It can proximate intergral of g(x) = g_tilde(x) * exp(-x^2) from -inf to inf

# get weights and nodes
x, w = gausshermite(3)

# suppose we want to calculate the expectation of f(x) where x is standard normally distributed: ∫f(x)*1/(sqrt(2π))*exp(-x^2/2)dx
# we need to do a change of variables

g_tilde(x,f) = f(sqrt(2)x)/sqrt(π)

# approximate the expectation with Gauss-Hermite Quadrature
E = sum(w.*g_tilde.(x,f))

# compare with Monte Carlo simulation approximation
val = 0.0
sims = 10^6
for i = 1:sims
    val += f(rand(Normal())) / sims
end
val