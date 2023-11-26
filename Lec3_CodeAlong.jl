using Plots, Optim
# Optim gives optimize function for optimization algorithms

############# 1. Univariate Box Constrained optimization ##########

f(x,y) = (x-y)^2

# Brent's method (this is default for univariate box constrained)

x -> f(x,1.0) # this creates an inline function

# save results as opts
opt = optimize(x -> f(x,1.0),
                        -5.0, # lower bound
                         5.0) # upper bound
                         
############# 2. Multivariate Case Without Derivative ###########

# Rosenbrock function
function Rosenbrock(x::Vector{Float64})
    (1-x[1])^2 + 100.0*(x[2] - x[1]^2)^2
end

Rosenbrock([1.0; 1.0]) # global minimum

# evaluate at a bunch of points so we can plot 
x_grid = collect(-3.0:0.01:3.0)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i in eachindex(x_grid), j in eachindex(x_grid)
    z_grid[i,j] = Rosenbrock([x_grid[i]; x_grid[j]])
end

# plot the Rosenbrock function
surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50,50))
contourf(x_grid, x_grid, log.(1.0 .+ z_grid), seriescolor=:inferno)

# Nelder-Mead
guess = [10.0; 10.0]
opt = optimize(Rosenbrock, guess)
opt.minimizer # this is the minimizer
opt.minimum

############# 2. Multivariate Case With Derivatives ###########

## gradient of Rosenbrock function
function gradient(G, x::Vector{Float64})
    G[1] = -2.0*(1.0 - x[1]) - 400.0*x[1]*(x[2] - x[1]^2)
    G[2] = 200.0*(x[2] - x[1]^2)
    return G
end

## Hessian of Rosenbrock function
function hessian(H, x::Vector{Float64})
    H[1,1] = 2 - 400.0*x[2] + 1200.0*x[1]^2
    H[1,2] = -400.0 * x[1]
    H[2,1] = -400.0 * x[2]
    H[2,2] = 200.0
    return H
end

# Newton's method is default when providing gradient and Hessian
guess = [0.0; 0.0]
opt = optimize(Rosenbrock, gradient, hessian, guess)
opt.minimizer
opt.minimum

# many minima case 

function Greiwank(x::Vector{Float64})
    val = (1/4000)*sum(x.^2) - prod(cos.(x./sqrt(length(x)))) + 1
    return val
end

## evaluate at a bunch of pts and plot 
x_grid = collect(-5.0:0.01:5)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i in 1:nx, j = 1:nx
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Greiwank(guess)
end

## plots 
surface(x_grid, x_grid, z_grid, seriescolor=:viridis, camera = (50, 70))
contourf(x_grid, x_grid, z_grid, seriescolor=:inferno)

### global optimum at (0,0)
guess_init = [3.0, 3.0]
opt = optimize(Greiwank, guess_init)
opt.minimizer
opt.minimum

guess_init = [1.0, 1.0]
opt = optimize(Greiwank, guess_init)
opt.minimizer
opt.minimum

# we can try multiple guesses to starting
function MultiStart()
    guess_grid = collect(-5:2.0:5)
    nx = length(guess_grid)
    minimum, minimizers = 100, [100, 100]

    for i = 1:nx, j = 1:nx
        guess = [guess_grid[i],guess_grid[j]]
        opt = optimize(Greiwank, guess)
        if opt.minimum < minimum
            minimum = opt.minimum
            minimizers = opt.minimizer
        end
    end

    return minimum, minimizers

end

min, minimizers = MultiStart()

############# 3. OLS Example ###########

using Distributions, Random

dist = Normal(0,1)
β_0 = 1.0
β_1 = 2.0
β_2 = 3.0
n = 10000
x = rand(n).*10
x2 = x.^2
Random.seed!(1234)
ϵ = rand(dist, n)
Y_true = β_0 .+ β_1.*x + β_2.*x2 .+ ϵ
X = hcat(ones(n), x, x2)
β_ols = inv(X'*X)*X'*Y_true

# estimate β using algorithms like Nelder-Mead
function sq_error(β::Vector{Float64})
    β_0, β_1, β_2 = β[1], β[2], β[3]
    Random.seed!(1234) # if we don't use it, it will fail
    ϵ = rand(Normal(), n)
    Y_true = 1.0 .+ 2.0.*x + 3.0.*x2 + ϵ
    Y_predict = β_0 .+ β_1.*x + β_2.*x2
    error = sum((Y_true.-Y_predict).^2)
    return error
end

# do OLS with Nelder-Mead
guess_init = [0.0, 0.0, 0.0]
opt = optimize(sq_error, guess_init)
opt.minimizer






