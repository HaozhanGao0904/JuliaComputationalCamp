
################################################################################
# Packages
################################################################################

using Optim, Interpolations, Plots



################################################################################
# Question 1
################################################################################

function Himmelblau(guess)
    x, y = guess[1], guess[2]
    val = (x^2 + y - 11)^2 + (x + y^2 - 7)^2 #compute himmelbau function
    val
end

#define x-y grids and compute function at each point
x_grid, y_grid = collect(-5.0:0.01:5.0), collect(-5.0:0.01:5.0)
nx, ny = length(x_grid), length(y_grid)
z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    x, y = x_grid[i], y_grid[j]
    z_grid[i,j] = log(Himmelblau([x, y])+1)
end

contourf(x_grid, y_grid, z_grid, seriescolor=:viridis) #plot
#savefig("ps2_3a.png")

surface(x_grid, y_grid, z_grid, seriescolor=:cinferno, camera = (50,30)) #plot
#savefig("ps2_3b.png")

#gradient
function g!(G, guess)
    x, y = guess[1], guess[2]
    G[1] = 2*(x^2 + y - 11) * 2*x + 2 * (x + y^2 - 7)
    G[2] = 2 * (x^2 + y - 11) + 2 * (x + y^2 - 7) * 2*y
end

#Hessian
function h!(H, guess)
    x, y = guess[1], guess[2]
    H[1] = 12*x^2 + 4*y - 44 + 2
    H[2] = 4*x + 4*y
    H[3] = 4*x + 4*y
    H[4] = 2 + 4*x + 12*y^2 - 28
end

#experimenting with algorithms and starting guesses
x_0 = [-100.0, 100.0]
@elapsed opt = optimize(Himmelblau, g!, h!, x_0)
@elapsed opt = optimize(Himmelblau, x_0)



################################################################################
# Question 2
################################################################################

function Ackley(guess)
    x, y = guess[1], guess[2]
    val = -20 * exp(-0.2 * sqrt(0.5 * (x^2 + y^2))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + ℯ + 20 #compute ackley
end

#define x-y grids and compute function at each point
x_grid, y_grid = collect(-4.0:0.01:4.0), collect(-4.0:0.01:4.0)
nx, ny = length(x_grid), length(y_grid)
z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    x, y = x_grid[i], y_grid[j]
    z_grid[i,j] = Ackley([x,y])
end

surface(x_grid, y_grid, z_grid, seriescolor=:darkrainbow, camera = (50,50)) #plot
#savefig("ps2_4a.png")

contourf(x_grid, y_grid, z_grid, seriescolor=:kdc) #plot
#savefig("ps2_4b.png")


#experimenting with algorithms and starting guesses
x_0 = [10.0,10.0]
opt = optimize(Ackley, x_0, LBFGS())
opt = optimize(Ackley, x_0)



################################################################################
# Question 3
################################################################################

function Rastrigin(guess::Vector{Float64})
    n = length(guess)
    val = 10 * n
    for i = 1:n
        val += guess[i]^2 - 10 * cos(2 * pi * guess[i]) #compute rastrigin function
    end
    val
end

#define x-y grids and compute function at each point
x_grid = collect(-5.12:0.01:5.12)
nx = length(x_grid)
rast_grid = zeros(nx)
for i = 1:nx
    guess = [x_grid[i]]
    rast_grid[i] = Rastrigin(guess)
end

plot(x_grid, rast_grid)
#savefig("ps2_5_first.png")


#2 dimensional
x1_grid, x2_grid = collect(-5.12:0.01:5.12), collect(-5.12:0.01:5.12)
nx = length(x1_grid)
rast_grid = zeros(nx, nx)
for i = 1:nx, j = 1:nx
    guess = [x1_grid[i], x2_grid[j]]
    rast_grid[i,j] = Rastrigin(guess)
end

surface(x1_grid, x2_grid, rast_grid, seriescolor=:fire, camera = (50,50)) #plot
#savefig("$dir/ps2_5a.png")

contourf(x1_grid, x2_grid, rast_grid, seriescolor=:coolwarm) #plot
#savefig("$dir/ps2_5b.png")


#experimenting with algorithms and starting guesses
x_0 = [0.0, 0.0]
opt = optimize(Rastrigin, x_0, LBFGS())
opt = optimize(Rastrigin, x_0)



################################################################################
# Question 4
################################################################################

function lin_approx(f, a::Float64, b::Float64, n::Int64, x::Float64)
    grid = collect(range(a, length = n, stop = b))
    func_grid = zeros(n)

    #fill in function grid
    for i = 1:n
        func_grid[i] = f(grid[i])
    end

    #point in discretized domain that is bigger than x
    ub_index = findfirst(z->z>x, grid)
    lb_index = ub_index - 1

    #compute interpolation
    interp = func_grid[lb_index] + (x - grid[lb_index]) * (func_grid[ub_index] - func_grid[lb_index])/(grid[ub_index] - grid[lb_index])
    interp #return
end

f(x) = x^2
lin_approx(f, 0.0, 5.0, 6, 2.5)



################################################################################
# Question 5
################################################################################

using Interpolations, Plots, Optim

function approx_log(x_grid::Vector{Float64})
    
    f(x) = log(1+x)

    f_grid = f.(x_grid)

    f_interp = LinearInterpolation(x_grid, f_grid)

    f_approx = f_interp.(collect(0:1000)./10)

    f_true = f.(collect(0:1000)./10)

    errors = abs.(f_approx - f_true)

    errors, f_approx, f_true #return
end

#initialize grids and compute approximation error
x_grid = collect(0.0:10.0:100.0)

x_fine = collect(0:0.1:100.0)

errors, f_approx, f_true = approx_log(x_grid)

sum(abs.(errors)) #report total error

plot(x_fine, errors) #plot error
#savefig("ps2_2a.png")
plot(x_fine, [f_approx, f_true]) #plot
#savefig("ps2_2b.png")

#optimize point placement
function eval_points(diffs::Vector{Float64})

    #Rule out negative step sizes
    if sum(diffs.<=0.0) > 0
        return Inf
    else

        grid = [0.0; [sum(diffs[1:i]) for i in 1:9]; 100.0]
        
        #rule out steps size too big
        if maximum(grid)>100
            return Inf
        else
            errs, func_approx, func_true = approx_log(grid) #compute approximation error
            return sum(errs) #return
        end

    end
end

#run optimizing algorithm
diffs_init = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
opt = optimize(diffs->eval_points(diffs), diffs_init; g_tol = 1e-4) #optimize

#report new results
grid_ans = [0.0; [sum(opt.minimizer[1:i]) for i in 1:9]; 100.0]
errs, func_approx, func_true = approx_log(grid_ans)
plot(x_fine, errs) #plot error
#savefig("$dir/ps2_2c.png")
plot(x_fine, [func_approx, func_true])
#savefig("$dir/ps2_2d.png")



