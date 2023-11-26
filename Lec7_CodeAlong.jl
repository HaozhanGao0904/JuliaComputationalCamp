################################################################################
# Set-Up for Multi-Processing
################################################################################

# Package for multi-processsing functions and marcros:
using Distributed

#Add 2 processes/workers
addprocs(2)



################################################################################
# Memory across processes
################################################################################

#Define r on master processes
r = rand()

#Other workers don't know what r is. 
@everywhere @show r #Error

#Define a different r on every process
@everywhere r = rand()
@everywhere @show r #Every process generated a different random number!



################################################################################
# Simple Example
################################################################################

#Suppose we have a vector of numbers x = [x₁, x₂, ..., xₙ]
#We want to calculate the sum of their square roots S = Σ(xᵢ^0.5)

@everywhere function sqrt_sum(A)
    S = 0.0
    for i in eachindex(A)
        #sleep(0.001)
        S += sqrt(A[i])
    end
    return S
end

#We can run this just on our master process
A = rand(1000)
@time sqrt_sum(A)

#We can also break the sum up into batches and send to our workers using @distributed. 
#The @sync macro tells the master to wait until the workers are done with the for loop. 
function sqrt_sum_distributed(A, N_batches)

    N = length(A)
    batch_size = Int(length(A) / N_batches)

    # The (+) in front of the for loop will add together the results of the calculations
    # inside the for loop and saves it as S
    S = @sync @distributed (+) for batch in [(1:batch_size) .+ offset for offset in 0:batch_size:(N-1)]
                            sqrt_sum(A[batch]) # Returns sum of this batch which will get added to S
    end
    
    return S

end

#This creates a vector of ranges
N = length(A)
batch_size = 100
[(1:batch_size) .+ offset for offset in 0:batch_size:(N-1)]

#Run distributed
@time sqrt_sum_distributed(A, 100)



#We can also run this in parellel using pmap()
function sqrt_sum_pmap(A, N_batches)

    N = length(A)
    batch_size = Int(length(A) / N_batches)

    S = sum(
            pmap(
                batch -> sqrt_sum(A[batch]), #Operation to be done by workers 
                [(1:batch_size) .+ offset for offset in 0:batch_size:(N-1)] #Vector of inputs to be broken up accross workers
                )
            )

    return S
end

@time sqrt_sum_pmap(A, 100)


@sync @distributed for i = [6; 2; 2; 2]
    @show i, sleep(i)
end

pmap(i -> (sleep(i); println(i)), [6; 2; 2; 2])



#=
Here, because the calculations we are doing inside the for loop (sqrt and add), 
there is not much speed gain from parallelization. Although we can go through these
calculations twice as fast on our two workers, there is an "overhead" cost of 
passing information back and forth between master and worker. Paralleization has
performance advantage when this overhead cost is cost is small relative to the computation
time inside the loop. 

If you uncomment the sleep command inside sqrt_sum(), you will see that @distributed and pmap
are now about twice as fast as just runnning on the master process. 

pmap is better when each iteration of the for loop takes a longer time. This is because pmap does
dynamic scheduling, and each worker will get a new task when they finish their last task. 

@distributed tends to run faster than pmap when the calculations inside the for loop are quick. This
is because @distributed does not do dynamic scheduling and just sends an equal number of tasks to each
processor at the start. Dynamic scheduling takes time to do, which is only worthwhile if the computations
completed in the for loop take a long time to run. 

=#

# I have two worker processes right now. 
#Each worker gets two tasks. This takes about 4 seconds
@time @sync @distributed for i = [3; 1; 1; 1]
    @show i, sleep(i)
end

# One worker gets the 3 second task. The other worker gets the two second tasks. 
# This takes only 3 seconds
@time pmap(i -> (sleep(i); println(i)), [3; 1; 1; 1])



####################################################################################################
# Dynamic programming examples
####################################################################################################



####################################################################################################
# Fibonacci numbers
####################################################################################################


# Naive Implementation.
function fib_naive(n::Int64)
    if n <= 1
        return n
    else
        return fib_naive(n-1) + fib_naive(n-2)
    end
end



# fib_naive(4)
# fib_naive(3) + fib_naive(2)
# fib_naive(2) + fib_naive(1) + fib_naive(1) + fib_naive(0)
# fib_naive(1) + fib_naive(0) + fib_naive(1) + fib_naive(1) + fib_naive(0) = 3

# fib_naive(6)
# fib_naive(5) + fib_naive(4)
# fib_naive(4) + fib_naive(3) + fib_naive(4) uh oh! we will repeat operations here by calculating 
#                                            the fourth fibonacci number twice



#counting up
function fib_dp(n::Int64)

    fib_vec = zeros(n+1)
    fib_vec[2] = 1

    for i = 3:n+1
    
        fib_vec[i] = fib_vec[i-1] + fib_vec[i-2]
    
    end

    return fib_vec[n+1]

end




#Fast:
@time fib_dp(42)
#Slow
@time fib_naive(42)

@time fib_dp(1000)
#@time fib_naive(100) #don't even try to run this.



####################################################################################################
# Egg Dropping
####################################################################################################

# Suppose we have n eggs and k floors. Where should we drop egg from to determine maximum drop of 
# egg fastest?


#Naive
function egg_drop_naive(n::Int64, k::Int64)

    # If we have one egg or there is only one floor. 
    if k == 1 || k == 0 || n == 1
        
        return k
    
    else
        
        candidate_min = 1e10 #something bad
        
        for x = 1:k #loop over floors to drop from
        
            val = max(egg_drop_naive(n-1, x-1), egg_drop_naive(n, k-x))
        
            if val<candidate_min
        
                candidate_min = val
        
            end
        
        end
        
        return candidate_min + 1
    
    end

end

# Very very very slow:
# egg_drop_naive(2,36)

####Dynamic Programming
function egg_drop_dp(n::Int64, k::Int64)
    
    trials = zeros(n, k) #matrix that stores number of trials required for any combination of eggs/floors up through and including [n,k]

    trials[:,1] .= 1

    for i = 1:k
        trials[1,i] = i
    end

    for  k_c = 2:k #Current number of floors left
    for  n_c = 2:n #Current number of eggs left.

        candidate_min = 1e10

        for j = 1:k_c

            val_1 = 0 # Number of additional trials if egg breaks. If floor = 1, then done.
            if j>1
                val_1 = trials[n_c-1,j-1]
            end

            val_2 = 0 #number of additional trials if egg does not break. If j = k_c, then done.
            if j < k_c
                val_2 = trials[n_c, k_c - j]
            end

            # Add one because we drop an egg this round
            res = 1 + max(val_1, val_2)
            if res<candidate_min
                candidate_min = res
            end
        
        end
    
        trials[n_c, k_c] = candidate_min
    
    end
    end
    
    trials

end

#Much faster
egg_drop_dp(4, 8)

####How???
#First trial: drop from floor 8. If egg breaks, then 1 egg and 7 flooors
#If not break: drop from floor 15. If breaks, check floors 9-14 with one egg.
#If not break: floor 21
#If not break: floor 26
#If not break: floor 30
#If not break: floor 33
#If not break: floor 35
#If not break: floor 36



####################################################################################################
# Minimum Travel Time between Two Cities
####################################################################################################

# What is the minimum travel time between cities two cities?
# The symmetric matrix d gives the distance between N cities. 
# We are interested in the distance from 1 to the other cities. 

d = [0.0 4.0 1.0; 
     4.0  0.0 1.0;
     1.0  1.0 0.0]

function min_distance(d)

    #Number of cities
    N = size(d)[1]

    min_dist = ones(N,N).* Inf
    visited = []

    #Visit source city
    min_dist[1,:] = d[1,:]
    visited = [1]

    for i in 2:N

        # List of unvisited cities
        unvisited = setdiff(1:N, visited)
        
        # Visit the city that has the least distance from source
        v = unvisited[argmin(min_dist[i-1,unvisited])]



        for j in 1:N
            min_dist[i,j] = minimum([min_dist[i-1,j];   # Go to city j without stopping in city v on the way
                                  min_dist[i-1,v] + d[v,j]] # Go to city j by going to city v and then going to city j. 
                                  )
        end

        visited = [visited; v]
    end

    min_dist

end


d = [0 4 2 1; 
        4 0 1 2; 
        2 1 0 0.5; 
        1 2 0.5 0
        ]
        
min_distance(d)



