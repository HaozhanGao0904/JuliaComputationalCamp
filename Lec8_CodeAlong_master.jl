################################################################################
# Packages
################################################################################

using Distributions
using Interpolations
using Optim
using Parameters
using Plots
using Random



################################################################################
# Load Structs and Functions
################################################################################

include("CodeAlong8_structs.jl")
include("CodeAlong8_functions.jl")



################################################################################
# Run Code
################################################################################

prim, res = solve_model()

I_avgs = zeros(51)
for i = 1:prim.N_t
    I_avgs[i] = mean(res.I[:,:,i])
end

plot(prim.t_grid, I_avgs)

h_avg, earnings_avg, I_avg, h_std, earnings_std, I_std  = Simulate(prim, res, 110.0, 7, 10000)

plot(prim.t_grid,h_avg,grid=false,ribbon=h_std,fillalpha=.25,linewidth=:4, ylim = [0,300])
plot(prim.t_grid,earnings_avg,grid=false,ribbon=earnings_std,fillalpha=.25,linewidth=:4, ylim = [0,250])
plot(prim.t_grid,I_avg,grid=false,ribbon=I_std,fillalpha=.25,linewidth=:4)

    
