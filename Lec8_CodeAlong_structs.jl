

################################################################################
# Model Primitives
################################################################################

@with_kw struct Primitives
   
    β::Float64 = 0.95 #discount rate
    ρ::Float64 = 2.0 #CRRA coefficient
    κ::Float64 = 0.7 #Ben-Porath HC coefficient

    #distribution of HC shocks
    μ::Float64 = -0.029
    σ::Float64 = 0.111
    N_ϵ::Int64 = 10
    ϵ_grid::Vector{Float64} = discretize_lognormal(μ, σ, N_ϵ)
    
    #grids and dimensions
    a_grid::Vector{Float64} = collect(range(0.01, length = 10, stop = 0.6)) #ability grid
    N_a = length(a_grid) #number of ability points

    h_grid::Vector{Float64} = collect(range(50, length = 50, stop = 400)) #HC grid
    N_h = length(h_grid)
    
    t_grid = collect(20:1:70)
    N_t = length(t_grid)

end

function discretize_lognormal(μ::Float64, σ::Float64, n::Int64)
    
    dist = LogNormal(μ, σ) #set up log normal distribution
    shocks = zeros(n) #preallocate market luck states

    for i = 1:n #begin loop to fill discretized vector
        quant = (2*i-1)/(2*n) #desired quantile
        shocks[i] = quantile(dist, quant) #fill element of state vector
    end

    shocks #return deliverable

end



################################################################################
# Model Results
################################################################################

@with_kw struct Results

    V::Array{Float64,3} #Value function
    I::Array{Float64,3} #Investment decision

end


