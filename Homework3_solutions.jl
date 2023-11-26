################################################################################
#Packages
################################################################################

using CSV
using DataFrames 
using Distributions
using GLM #Regressions
using Optim
using Plots
using Random
using StatsBase



################################################################################
# Question 1
################################################################################

#Set working directory where data is
cd("C:/Users/kghun/Downloads")

#Read in dataset 
df = CSV.read("lwage.csv", DataFrame)

#Create quadratic exp
df.exp_sq = df.exp .^2

#Run OLS
lmsum = lm(@formula(lwage ~ educ + exp + exp_sq), df)

coef_OLS = coef(lmsum)
σ_OLS = std(residuals(lmsum))



#Log Likelihood function
function loglikelihood_LR(n, X, Y, θ)

    #n = length(Y)

    #Unpack θ
    #β_0 = θ[1]
    #β_1 = θ[2]
    #β_2 = θ[3]
    #β_3 = θ[4]
    #σ = θ[5]

    #Log Likelihood
    -n/2*log(2π) - n/2* log(θ[5]^2) - (sum((Y .- θ[1] - X * θ[2:4] ).^2) / (2*θ[5]^2)) 

end


# Optimize Log Likelihood with initial guess as OLS result
Y = Vector{Float64}(df.lwage)
X = Matrix{Float64}(df[:,2:4])
opt = optimize(θ -> -loglikelihood_LR(length(Y),X, Y, θ), [coef_OLS; σ_OLS])

#MLE results are same as OLS
coef_MLE = opt.minimizer[1:4]
σ_MLE = opt.minimizer[5]





################################################################################
# Question 2
################################################################################

function get_matches(n::Int64, sims::Int64)
    
    matches = zeros(sims)
    
    for i = 1:sims
        
        people = collect(1:1:n)
        slips = collect(1:1:n)
        
        slips = shuffle(slips) #shuffle slips
        
        matches[i] = sum(people .== slips) #update MC results vector
    
    end
    
    matches

end

#obtain results and plot
matches10 = get_matches(10, 10000)
matches20 = get_matches(20, 10000)

histogram(matches10)
histogram(matches20)




################################################################################
# Question 3
################################################################################


function Retirement(save_rate, sims, risk; 
                    μ_r = 0.06, σ_r = 0.06, raise_lb = 0.0, raise_ub = 0.06, 
                    s = 20230813)

    Random.seed!(s)
    
    dist_r = Normal(μ_r, σ_r) #distribution of investment return
    dist_raise = Uniform(raise_lb, raise_ub) #Distribution of raises
    μ_raise = (raise_lb + raise_ub ) / 2 #Average Raise

    results = zeros(sims)

    for i = 1:sims

        #earnings/savings at 30
        earnings = 100
        savings = 100

        #Loop over all ages
        for j = 30:67

            raise = ifelse(risk, rand(dist_raise), μ_raise)
            r = ifelse(risk, rand(dist_r), μ_r)

            savings *= (1 + r) #adjust savings
            earnings *= (1 + raise) #adjust earnings
            
            savings += (earnings * save_rate) #add to savings

        end

        results[i] = savings/earnings #divide by salary
    
    end
    
    results

end

#Evaluate once
res_norisk = Retirement(.1125, 1, false)
mean(res_norisk)

#Find savings save_rate that gives 10x earnings with no risk
opt = optimize(savings_rate -> (Retirement(savings_rate, 1, false)[1] - 10.0)^2, 0.0, 20.0)
opt.minimizer
Retirement(opt.minimizer, 1, false)[1]

# Simulate with risk
res_risk = Retirement(opt.minimizer, 10000, true)
histogram(res_risk)
mean(res_risk .< 10)

# Find savings rate that gives 90% of getting 10x earnings
opt = optimize(savings_rate -> (mean(Retirement(savings_rate, 10000, true).>=10.0)-0.9)^2, 0.0, 0.2)
opt.minimizer
res = Retirement(opt.minimizer, 10000, true)
histogram(res)



