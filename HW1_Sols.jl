############### Q1 ###############
function factorial1(n::Int64)
    value = 1
    while n>0
        value *= n
        n -= 1
    end
    return value
end

factorial1(5)

function factorial2(n::Int64)
    value = 1
    for i = 1:n
        value *= i
    end

    return value
end

factorial2(5)

############# Q2 #################

function p(x::Float64, coeff::Vector{Float64})
    val = 0
    for (i, coef) in enumerate(coeff)
        val += coef * x^(i-1)
    end
    return val
end

p(1.0, [1.0, 2.0, 3.0, 4.0]) # must match the type

############# Q3 ##############

function approx_π(n::Int64)
    x = rand(n)
    y = rand(n)
    num = 0
    for i in 1:n
        if x[i]^2 + y[i]^2 <= 1
            num += 1
        end
    end
    π = 4 * num / n
end

approx_π(10000000)

############### Q4 ################

using Distributions, Plots, Parameters
function ols(n::Int64, m::Int64)
    a, b, c, d, σ = 0.1, 0.2, 0.5, 1.0, 0.1
    y = zeros(n)
    coeff_vec = zeros(4, m)
    x_1, x_2 = rand(Normal(), n), rand(Normal(), n)
    for i in 1:m
        w = rand(Normal(), n)
        X = hcat(x_1, x_1.^2, x_2, ones(n))
        y = a.*x_1 .+ b.*x_1.^2 .+ c.*x_2 .+ 1.0 .+ σ.*w
        coeff = inv(X'*X)*X'*y
        coeff_vec[:, i] = coeff
    end
    return coeff_vec
end

ols_vec = ols(1000, 1000)

plot_a = histogram(ols_vec[1,:], title = "Estimates of a")
plot_b = histogram(ols_vec[2,:], title = "Estimates of b")
plot_c = histogram(ols_vec[3,:], title = "Estimates of c")
plot_d = histogram(ols_vec[4,:], title = "Estimates of d")

############## Q5 ###############

using Distributions, Plots

function simulate_walk(n::Int64, α::Float64, σ::Float64, t_max::Int64)
    first_times_to_cross = zeros(n)

    for i = 1:n # repeat simulation n times
        cross_yet = false
        x_now = 1
        for t = 1:t_max
            x_next = α*x_now + σ*rand(Normal())
            if x_next <= 0
                first_times_to_cross[i] = t
                cross_yet = true
                break
            end
            x_now = x_next
        end
        
        if !cross_yet
            first_times_to_cross[i] = t_max
        end
    end
    
    return first_times_to_cross
end

# α = 0.8
temp = simulate_walk(100, 0.8, 0.2, 200)
histogram(temp)
# α = 1.0
temp = simulate_walk(100, 1.0, 0.2, 200)
histogram(temp)
# α = 1.2
temp = simulate_walk(100, 1.2, 0.2, 200)
histogram(temp)

############## Q6 ##################

function Newton(f, fprime, x_0::Float64, tol::Float64, maxiter::Int64)\
    error,root = 100, 0 # define initial error and root
    iter = 1
    while error > tol
        root = x_0 - f(x_0)/fprime(x_0)
        error = abs(root-x_0)
        x_0 = root
        iter += 1
        if iter > maxiter
            print("have reached the maximum iteration \n")
            break
        end
    end
    return root 
end

# take f(x) = (x-1)^3 as example
f(x) = (x-1)^3
fprime(x) = 3(x-1)^2
root = Newton(f, fprime, 4.0, 1e-3, 100)

# different function: f(x) = ln(x) + 3x - 7
f(x) = log(x) + 3x - 7
fprime(x) = 1/x + 3
root = Newton(f, fprime, 4.0, 1e-3, 100)

########### Q7 ############
using Parameters, Plots

# struct for model parameters
@with_kw struct ModelParameters
    β::Float64 = 0.99 
    δ::Float64 = 0.025
    α::Float64 = 0.36
    k_grid::Vector{Float64}=collect(range(0.1, length = 1000, stop = 45.0))       
    N_k::Int64 = length(k_grid)
    z_grid::Vector{Float64} = [1.25; 0.2]
    N_z::Int64 = 2
    Π::Array{Float64,2} = [0.977 0.023; 0.074 0.926]
    tol::Float64 = 10^-4
end

# struct for model solutions
@with_kw struct ModelSolutions
    V::Array{Float64,2}
    kp::Array{Float64,2}
end

function build_ModelSolutions(params)
    V = zeros(params.N_k, params.N_z)
    kp = zeros(params.N_k, params.N_z)

    sols = ModelSolutions(V, kp)

    return sols
end

function build_structs()
    params = ModelParameters()
    sols = build_ModelSolutions(params)

    return params, sols
end

# Bellman operator
function Bellman(params, sols)
    @unpack_ModelParameters params
    @unpack_ModelSolutions sols

    V_next = zeros(N_k,N_z)
    kp_next = zeros(N_k,N_z)

    for i_k = eachindex(k_grid), i_z = eachindex(z_grid)
        max_u = -1e10
        k = k_grid[i_k]
        z = z_grid[i_z]
        budget = z*k^α + (1-δ)*k

        for i_kp = eachindex(k_grid)

            c = budget - k_grid[i_kp]
            
            if c > 0
                V_temp = log(c) + β*(Π[i_z,1]*V[i_kp,1] + Π[i_z,2]*V[i_kp,2])

                if V_temp > max_u
                    max_u = V_temp
                    kp_next[i_k,i_z] = k_grid[i_kp]
                end

            end
        end

        V_next[i_k,i_z] = max_u

    end

    return V_next, kp_next

end

# solve model
function solve_model(params, sols)
    @unpack_ModelParameters params
    @unpack_ModelSolutions sols

    V_next = zeros(N_k,N_z)
    kp_next = zeros(N_k,N_z)
    max_diff = tol + 10.0
    n = 0

    while max_diff > tol
        n += 1
        V_next, kp_next = Bellman(params, sols)

        max_diff = maximum(abs.(V_next-V))
        V .= V_next
        kp .= kp_next

        @show n, max_diff
    end
end

params, sols = build_structs()

@elapsed solve_model(params, sols)

plot(params.k_grid, sols.V)
plot(params.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))


                
            


