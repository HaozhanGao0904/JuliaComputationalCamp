using Optim, Interpolations, Plots

############ Q1 ############

function Himmelblau(guess)
    x,y = guess[1], guess[2]
    val = (x^2 + y - 11)^2 + (x + y^2 -7)^2
    return val
end

# define x-y grids
x_grid, y_grid = collect(-5.0:0.01:5.0), collect(-5.0:0.01:5.0)
nx, ny = length(x_grid), length(y_grid)
z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    x, y = x_grid[i], y_grid[j]
    z_grid[i,j] = log(Himmelblau([x,y]))
end

contourf(x_grid, y_grid, z_grid, seriescolor=:viridis)

surface(x_grid, y_grid, z_grid, seriescolor=:inferno, camera=(50,20))

# gradient
function gradient(G, guess)
    x,y = guess[1], guess[2]
    G[1] = 2*(x^2 + y - 11) * 2*x + 2 * (x + y^2 - 7)
    G[2] = 2 * (x^2 + y - 11) + 2 * (x + y^2 - 7) * 2*y
end

# Hessian
function Hessian(H, guess)
    x, y = guess[1], guess[2]
    H[1] = 12*x^2 + 4*y - 44 + 2
    H[2] = 4*x + 4*y
    H[3] = 4*x + 4*y
    H[4] = 2 + 4*x + 12*y^2 - 28
end

x_0 = [-100.0, 100.0]

# Newton's method
@elapsed opt = optimize(Himmelblau, gradient, Hessian, x_0)
opt.minimizer
# Nelder-Mead
@elapsed opt = optimize(Himmelblau,x_0)
opt.minimizer

############ Q2 ############
function Ackley(guess)
    x,y = guess[1], guess[2]
    val = -20 * exp(-0.2 * sqrt(0.5 * (x^2 + y^2))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + ℯ + 20
    # type \euler for constant ℯ
end

x_grid, y_grid = collect(-4.0:0.01:4.0), collect(-4.0:0.01:4.0)
nx, ny = length(x_grid), length(y_grid)
z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    x, y = x_grid[i], y_grid[j]
    z_grid[i,j] = Ackley([x,y])
end

surface(x_grid,y_grid,z_grid,seriescolor=:rainbow, camera=(50,30))
contourf(x_grid,y_grid,z_grid,seriescolor=:kdc)

x_0 = [10.0,10.0]
# using LBFGS
opt = optimize(Ackley, x_0, LBFGS())
opt.minimizer
# using Nelder-Mead
opt = optimize(Ackley, x_0)
opt.minimizer
# it seems that Nelder-Mead is better for finding global minimum

############# Q3 #############

function Rastrigin(guess)
    n = length(guess)
    val = 10*n
    for i = 1:n
        val += guess[i]^2 - 10*cos(2*π*guess[i])
    end
    return val
end

# part A
x_grid = collect(-5.12:0.01:5.12)
nx = length(x_grid)
f_grid = zeros(nx)

for i = 1:nx
    guess = [x_grid[i]]
    f_grid[i] = Rastrigin(guess)
end

plot(x_grid,f_grid)

# part B
x1_grid, x2_grid = collect(-5.12:0.01:5.12), collect(-5.12:0.01:5.12)
nx = length(x1_grid)
f_grid = zeros(nx, nx)

for i = 1:nx, j = 1:nx
    guess = [x1_grid[i],x2_grid[j]]
    f_grid[i,j] = Rastrigin(guess)
end

surface(x1_grid,x2_grid,f_grid,seriescolor=:fire,camera=(50,20))
contourf(x1_grid,x2_grid,f_grid,seriescolor=:coolwarm)

# Part C 
x_0 = [100.0,100.0]
# using LBFGS
opt = optimize(Rastrigin,x_0,LBFGS())
opt.minimizer

# using Nelder-Mead
opt = optimize(Rastrigin, x_0)
opt.minimizer

############ Q4 #############

function linear_approx(f, a::Float64, b::Float64, n::Int64, x::Float64)
    grid = collect(range(a, length=n, stop=b))
    fcn_grid = zeros(n)

    for i = 1:n
        fcn_grid[i] = f(grid[i])
    end

    # point in discretized domain that is bigger than x
    upperbdindex = findfirst(z->z>x,grid)
    lowerbdindex = upperbdindex-1

    interp = fcn_grid[lowerbdindex] + (x - grid[lowerbdindex])*(fcn_grid[upperbdindex]-fcn_grid[lowerbdindex])/(grid[upperbdindex]-grid[lowerbdindex])
    return interp
end

f(x) = x^3
linear_approx(f,0.0,5.0,6,3.0)

h(x) = x^2
linear_approx(h, 0.0, 5.0, 6, 2.5)

############ Q5 ############ 

# part A
function approx_log(x_grid::Vector{Float64})

    f(x) = log(1+x)

    f_grid = f.(x_grid)

    f_interp = LinearInterpolation(x_grid, f_grid)

    f_approx = f_interp.(collect(0:0.1:100))

    f_true = f.(collect(0:0.1:100))

    errors = abs.(f_approx - f_true)

    return errors, f_approx, f_true
end

x_grid = collect(0.0:10.0:100.0)

x_fine = collect(0.0:0.1:100.0)

errors, f_approx, f_true = approx_log(x_grid)

sum(abs.(errors))

plot(x_fine, errors)
# this shows errors occur more severely near 0
# since the derivative changes fast near 0

plot(x_fine, [f_approx, f_true])

# part B

function eval_points(diffs::Vector{Float64})

    # rule out negative step sizes
    if sum(diffs.<=0) > 0
        return Inf
    else
        grid = [0.0; [sum(diffs[1:i]) for i in 1:9]; 100.0]

        if maximum(grid)>100
            return Inf
        else
            errors, func_approx, func_true = approx_log(grid)
            return sum(errors)
        end
    end
end

# part C 

# optimizing grid
diffs_init = fill(10.0,9)
opt = optimize(diffs->eval_points(diffs), diffs_init; g_tol = 1e-4)
opt.minimizer

# new interpolation
grid_opt = [0.0; [sum(opt.minimizer[1:i]) for i in 1:9]; 100.0]
errors, func_approx, func_true = approx_log(grid_opt)
sum(errors)

# plot errors and interpolation function
plot(x_fine, errors)
plot(x_fine,[func_approx, func_true])














