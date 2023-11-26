###################################################################################################
#                                   Optimal Investment Problem                                    #
###################################################################################################

using Parameters 
#This package gives us the @with_kw macro, which allows us to define default values in our structs

### Struct for our model paramters
@with_kw struct ModelParameters

    β::Float64 = 0.99 #We can define default values because of the @with_kw macro. 
    δ::Float64 = 0.025
    α::Float64 = 0.36

    k_grid::Vector{Float64}=collect(range(0.1, length = 1800, stop = 45.0))       
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

    for i_k = eachindex(k_grid)
        
        max_util = -1e10
        k = k_grid[i_k]
        budget = k^α + (1-δ)*k

        for i_kp = eachindex(k_grid)
            
            c = budget - k_grid[i_kp]

            if c > 0
                
                V_temp = log(c) + β*V[i_kp]
            
                if V_temp > max_util
                    max_util = V_temp
                    kp_next[i_k] = k_grid[i_kp]
                end
                
            end

        end
        V_next[i_k] = max_util
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

@time sols = solve_model(para,sols)
 

using Plots
plot(para.k_grid, sols.V)
plot(para.k_grid, sols.kp)
plot!(collect(0:45), collect(0:45))


