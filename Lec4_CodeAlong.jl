# Packages
using Optim, Interpolations, Plots, Polynomials, Parameters

# polynomail approximation

# Runge Phenomenom

Runge(x) = 1/(1+25*x^2)
x_fine = collect(-1.0:0.01:1.0)
plot(x_fine, Runge.(x_fine))

# 5-degree polynomial approximation
x_coarse = collect(range(-1.0, 1.0, length=6))
poly_fit = fit(x_coarse, Runge.(x_coarse), 5)
plot!(x_fine, poly_fit.(x_fine)) # it didn't do a great job!
# the ! adds the plot/scatter on top of the previous plot rather than making a new plot

# let's try a higher degree polynomial
x_coarse = collect(range(-1.0, length = 10, stop = 1.0))
poly_fit = fit(x_coarse, Runge.(x_coarse), 9)
plot!(x_fine, poly_fit.(x_fine)) # this is even worse! there are waves on the tails


# interpolation

# linearly interpolation log function

x = 0.1:2:10.1
y = log.(x)

## make x_fine for plotting
x_fine = collect(0.1:0.1:10.1)
y_fine = log.(x_fine)

plot(x_fine, y_fine)
scatter!(x,y)

interp_y = LinearInterpolation(x, y)
plot!(x_fine, interp_y.(x_fine))

# one way to improve the interpolation is to have more points where the derivative is changing a LinearInterpolation
x_uneven = 0.1 .+ 10.0 .*collect(0.0:0.2:1.0).^2
y_uneven = log.(x_uneven)
interp_y_uneven = LinearInterpolation(x_uneven, y_uneven)
plot(x_fine, log.(x_fine))
plot!(x_fine, interp_y_uneven(x_fine))

# extrapolation: outside of your k_grid
interp_y(12) # this exeeds the grid
interp_y_extra = LinearInterpolation(x, y, extrapolation_bc = Line())
interp_y_extra(12)
log(12) # this is close since 12 is close to the grid

interp_y_extra(100)
log(100) # now this is not close since 100 is too far away

# cubic interpolation of the log function
x = 0.1:2:10.1
y = log.(x)
x_fine = collect(0.1:0.1:10.1)
y_fine = log.(x_fine)

cubic_y = CubicSplineInterpolation(x, y)
plot(x_fine, log.(x_fine))
plot!(x_fine, cubic_y.(x_fine))

# cubic spline is possible with uneven grids
# but Interpolation pacakge does not supprt uneven grid for cubic
# we can use Dierckx pacakge Instead

# extrapolation
cubic_y_extra = CubicSplineInterpolation(x, y, extrapolation_bc = Line())
plot!(x_fine, cubic_y_extra.(x_fine))

# bilinear interpolation
f(x,y) = 1 + x^2 + y^2
grid_coarse = collect(0.0:1.0:5.0)
grid_fine = collect(0.0:0.01:5.0)
n_coarse, n_fine = length(grid_coarse), length(grid_fine)
z_fine = zeros(n_fine, n_fine)

for i in eachindex(grid_fine), j in eachindex(grid_fine)
    z_fine[i,j] = f(grid_fine[i], grid_fine[j])
end

contourf(grid_fine, grid_fine, z_fine)

z_coarse = zeros(n_coarse, n_coarse)

for i in eachindex(grid_coarse), j in eachindex(grid_coarse)
    z_coarse[i,j] = f(grid_coarse[i], grid_coarse[j])
end

interp_z = LinearInterpolation((grid_coarse, grid_coarse), z_coarse)

grid_z_interp = zeros(n_fine, n_fine)
for i in eachindex(grid_fine), j in eachindex(grid_fine)
    grid_z_interp[i,j] = interp_z(grid_fine[i], grid_fine[j])
end
contourf(grid_fine, grid_fine, grid_z_interp)
contourf(grid_fine, grid_fine, z_fine .- grid_z_interp) # this is the difference between interpolation and the real 
contourf(grid_fine, grid_fine, (z_fine .- grid_z_interp)./z_fine) # this is the percentage difference between interpolation and the real 


# a better version of code of the best investment Problem
using Parameters
@with_kw struct ModelParameters

    β::Float64 = 0.99 #We can define default values because of the @with_kw macro. 
    δ::Float64 = 0.025
    α::Float64 = 0.36

    k_grid::Vector{Float64}=collect(range(0.1, length = 100, stop = 45.0)) #### Lower the Number of grid points to 100      
    N_k::Int64 = length(k_grid)

    tol::Float64 = 10^-4

end


@with_kw struct ModelSolutions

    V::Vector{Float64}
    kp::Vector{Float64}

end 

function build_ModelSolutions(para)

    V = zeros(para.N_k)
    kp = zeros(para.N_k)

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

    V_next = zeros(para.N_k)
    kp_next = zeros(N_k)

    #Interpolate value function for continuation value
    V_interp = LinearInterpolation(k_grid, V, extrapolation_bc = Line())

    for i_k = eachindex(k_grid)
        k = k_grid[i_k]
        budget = k^α + (1-δ)*k

        #Replace grid search with box constrained optimization!
        opt = optimize(kp -> -log(budget-kp) - β*V_interp(kp), 0.0, budget) # the last 2 numbers are the box constraint

        V_next[i_k] = -opt.minimum
        kp_next[i_k] = opt.minimizer

    end

    return V_next, kp_next

end



### Solve model
function solve_model(para, sols)
    para, sols = build_structs();    
    @unpack_ModelParameters para
    @unpack_ModelSolutions sols

    V_next = zeros(N_k)
    max_diff = tol + 10.0
    n = 0
    while max_diff > tol
        n +=1
        V_next, kp_next = bellman(para,sols)
        
        max_diff = maximum(abs.(V_next - V))
        sols.V .= V_next
        sols.kp .= kp_next

        @show n, max_diff

    end

    sols

end

para, sols = build_structs();

#So fast and just as accurate!
@time sols = solve_model(para,sols)
 

using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))