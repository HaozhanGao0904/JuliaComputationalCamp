################################################################################
# Simulation Function 
################################################################################

function Simulate(prim::Primitives, res::Results, h_0::Float64, i_a::Int64, N_sim::Int64)
    @unpack_Primitives prim
    @unpack_Results res

    Random.seed!(20230822)


    a = a_grid[i_a] #level of ability
    dist = LogNormal(μ, σ) #set up log normal distribution

    #simulate some outcomes for a representative person
    h_sims = zeros(N_t, N_sim) #tracking human capital
    earnings_sims = zeros(N_t, N_sim) #tracking earnings
    I_sims = zeros(N_t, N_sim) #tracking self-investment

    for s = 1:N_sim #loop over simulations
        
        h = h_0

        for i_t = 1:N_t

            I_interp = LinearInterpolation(h_grid, I[i_a,:,i_t], extrapolation_bc = Line())

            I_c = minimum([maximum([I_interp(h); 0.0]); 1.0])
            

            h_sims[i_t,s] = h
            earnings_sims[i_t, s] = h * (1-I_c)
            I_sims[i_t, s] = I_c

            hp = (h + a*(h * I_c)^κ) * rand(dist)
            h = hp

        end
    end

    h_avg, earnings_avg, I_avg = zeros(N_t), zeros(N_t), zeros(N_t)
    h_std, earnings_std, I_std = zeros(N_t), zeros(N_t), zeros(N_t)
    @show quantile(h_sims, [0.99])
    #a few lines to write here
    for i_t = 1:N_t
        h_avg[i_t] = mean(h_sims[i_t,:])
        earnings_avg[i_t] = mean(earnings_sims[i_t,:])
        I_avg[i_t] = mean(I_sims[i_t,:])
        h_std[i_t] = std(h_sims[i_t,:])
        earnings_std[i_t] = std(earnings_sims[i_t,:])
        I_std[i_t] = std(I_sims[i_t,:])
    end

    h_avg, earnings_avg, I_avg, h_std, earnings_std, I_std
end



################################################################################
#Function that solves the model
################################################################################

function solve_model()
    
    prim = Primitives()
    
    V = zeros(prim.N_a, prim.N_h, prim.N_t)
    I = zeros(prim.N_a, prim.N_h, prim.N_t)
    res = Results(V, I)
    
    backward_induct(prim, res)
    
    prim, res #return deliverables

end



################################################################################
#Backward induction protocol
################################################################################

function backward_induct(prim::Primitives, res::Results)
    
    @unpack_Primitives prim

    for i_t in reverse(eachindex(t_grid))
    
        println("Working on ", t_grid[i_t])
        res.V[:,:,i_t] .= bellman(prim, res, i_t)
    
    end
    
end



################################################################################
# Bellman Equation
################################################################################

function bellman(prim::Primitives, res::Results, i_t::Int64)
    @unpack_Primitives prim
    

    
    V_next = zeros(N_a, N_h) #value function to fill

    if i_t == N_t #check for maximal age
        
        for (i_h, h) in enumerate(h_grid)  #loop over state space
        
            V_next[:, i_h] .= u(h, ρ)
            res.I[:, i_h, i_t] .= 0.0 #no self-investment
        
        end
    
    elseif i_t != N_t #not in terminal age
    
        for (i_a, a) in enumerate(a_grid) 
            
            Vp_interp = LinearInterpolation(h_grid, res.V[i_a, :, i_t+1], extrapolation_bc = Line())
            
            for (i_h,h) in enumerate(h_grid) #loop over state space
        
                #Decide current investment I_c
                opt = optimize(I_c -> -objective(prim, Vp_interp, i_t, h, a, I_c), 
                                0.0, 1.0) 

                V_next[i_a, i_h] = -opt.minimum
                res.I[i_a, i_h, i_t] = opt.minimizer

            end
    
        end
    
    end
    
    V_next #return

end



################################################################################
# objective function
################################################################################

function objective(prim::Primitives, Vp_interp, i_t::Int64, h::Float64, a::Float64, I::Float64)
    @unpack_Primitives prim

    c = h * (1-I)
    val = u(c, ρ)

    # Take expectation. 
    for ϵ in ϵ_grid
    
        #Human capital next period hp
        hp = (h + a*(h*I)^κ) * ϵ
        val += (β * Vp_interp(hp))/N_ϵ
    
    end

    val

end



################################################################################
# CRRA Utility Function
################################################################################

function u(c::Float64, ρ::Float64)
    return (c^(1-ρ))/(1-ρ)
end


