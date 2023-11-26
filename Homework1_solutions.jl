################################################################################
# Problem 1
################################################################################

function factorial2(n::Int64)
    val = 1 #initialize value
    for i = 1:n
        val*=i #multiply by all integers between 1 and n
    end
    val #return
end



################################################################################
# Problem 2
################################################################################

function p(x::Float64, coeff::Vector{Float64})
    val = 0 #initializes
    for (i, coef) in enumerate(coeff)
        val += x^(i-1) * coef
    end
    val
end
p(1.0, [2.0, 3.0, 4.0])



################################################################################
# Problem 3
################################################################################

function approx_pi(iter::Int64)
    avg = 0
    for i = 1:iter
        x, y = rand(2) #bivariate uniform draw
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2) #distance from middle of square; pythagorean theorem
        avg+=(dist<0.5)
    end
    avg/=iter #approximation of area
    avg*=4 #π = A / r^2, where r = 1/2. So here π = 4A
    avg
end
approx = approx_pi(1000000)



################################################################################
# Problem 4
################################################################################

using Distributions, Plots

function simulate_ols(M::Int64, N::Int64)
    a, b, c, d, σ = 0.1, 0.2, 0.5, 1.0, 0.1 #true parameter values
    ols_vec = zeros(M, 4)
    x_1, x_2 = rand(Normal(), 50), rand(Normal(), 50) #random draws

    #loop over simulations
    for i = 1:M
        w = rand(Normal(), N) #randomly draw values of w

        #construct X and Y
        X = hcat(x_1, x_1.^2, x_2, ones(N))
        Y = a .* x_1 .+ b.*x_1.^2 .+ c.*x_2 .+ 1.0 .+ σ .* w
        ols = inv(X' * X) * X' * Y
        ols_vec[i,:] = ols
    end
    ols_vec
end

#Simulate and plot
ols_vec = simulate_ols(200,50)
p1 = histogram(ols_vec[:,1], title="Estimates of a")
p2 = histogram(ols_vec[:,2], title="Estimates of b")
p3 = histogram(ols_vec[:,3], title="Estimates of c")
p4 = histogram(ols_vec[:,4], title="Estimates of d")
plot(p1, p2, p3, p4, layout = (2,2), legend=:none)
#savefig("C://Users//kghun//OneDrive//Desktop//ComputationCamp//ps1_4.png")


################################################################################
# Problem 5
################################################################################

using Distributions, Plots

function simulate_walk(n::Int64, α::Float64, σ::Float64, t_max::Int64)
    first_times = zeros(n)

    for i = 1:n #loop over simulations
        cross_yet = false
        x_now = 1
        for t = 1:t_max
            x_next = α * x_now + σ * rand(Normal())
            if x_next<=0 #first crossing
                first_times[i] = t
                cross_yet = true
                break
            end
            x_now = x_next
        end
        if ! cross_yet #never crosses
            first_times[i] = t_max
        end
    end
    first_times #return
end

#0.8
temp = simulate_walk(100, 0.8, 0.2, 200)
histogram(temp)
#savefig("C:/Users/Garrett/Desktop/ps1_5a.png")

#1.0
temp = simulate_walk(100, 1.0, 0.2, 200)
histogram(temp)
#savefig("C:/Users/Garrett/Desktop/ps1_5b.png")

#1.2
temp = simulate_walk(100, 1.2, 0.2, 200)
histogram(temp)
#savefig("C:/Users/Garrett/Desktop/ps1_5c.png")



################################################################################
# Problem 6
################################################################################

function Newton(f, fp, x_0::Float64, tol::Float64, maxiter::Int64)
    error, root = 100, 0
    while error>tol
        root = x_0 - f(x_0)/fp(x_0)
        error = abs(root - x_0)
        x_0 = root
    end
    root
end

f(x) = (x-1)^3
fp(x) = 3 * (x-1)^2
root = Newton(f, fp, 4.0, 1e-3, 100)

#different function: f(x) = ln(x) + 3x - 7
f(x) = 2*log(x) + 3*x - 7
fp(x) = (2/x) + 3
root = Newton(f, fp, 4.0, 1e-3, 100)



###################################################################################################
# Question 7: Optimal Investment Problem
###################################################################################################

using Parameters 

### Struct for our model paramters
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

### Struct for our Model solutions
@with_kw struct ModelSolutions

    V::Array{Float64,2}
    kp::Array{Float64,2}

end 

function build_ModelSolutions(para)

    V = zeros(para.N_k, para.N_z)
    kp = zeros(para.N_k, para.N_z)

    sols = ModelSolutions(V,kp)

    return sols

end

function build_structs()

    para = ModelParameters()
    sols = build_ModelSolutions(para)

    return para, sols

end


### Bellman operator
function bellman(para, sols)
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    V_next = zeros(N_k,N_z)
    kp_next = zeros(N_k,N_z)

    for i_k = eachindex(k_grid), i_z = eachindex(z_grid)
        max_util = -1e10
         k = k_grid[i_k]
         z = z_grid[i_z]
        budget = z*k^α + (1-δ)*k

        for i_kp = eachindex(k_grid)
            
             c = budget - k_grid[i_kp]

            if c > 0
                
                 V_temp = log(c) + β*(Π[i_z,1]*V[i_kp,1] + Π[i_z,2]*V[i_kp,2])
            
                if V_temp > max_util
                    max_util = V_temp
                    kp_next[i_k,i_z] = k_grid[i_kp]
                end

            end

        end
         V_next[i_k,i_z] = max_util

    end

    return V_next, kp_next

end



### Solve model
function solve_model(para, sols)    
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    V_next = zeros(N_k,N_z)
    kp_next = zeros(N_k,N_z)
    max_diff = tol + 10.0
    n = 0

    while max_diff > tol
        n +=1
        V_next, kp_next = bellman(para,sols)
        
        max_diff = maximum(abs.(V_next - V))
        V .= V_next
        kp .= kp_next

        @show n, max_diff

    end

end

para, sols = build_structs();

@elapsed solve_model(para,sols)


using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))


