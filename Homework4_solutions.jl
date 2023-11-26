################################################################################
# Problem 1: Making Change
################################################################################

function get_changeways(S::Vector{Int64}, N::Int64)
    
    S = sort(S)
    size_S = length(S)
    
    table = zeros(size_S, N+1) # Columns are ways to make change for 0, 1, 2, ..., N.
    # The i,j entry of table is the number of ways to make change for j with the first i coins.

    #Only one way to make change for zero: Give no coins
    table[:,1] .= 1

    #Loop over all coins in set
    for (i, S_i) in enumerate(S)

        #Number of ways to make change without the current coin so far (previous row)
        if i > 1
            table[i,:] .= table[i-1,:]
        end
        
        #Loop over values up to N
        for val = 1:N

            #If value we want to make change for is at least as big as coin value:
            if val >= S_i
                
                #Add the number of ways we can make change for val - S_i to 
                #the number of ways we can make change for val. 
                table[i,val+1] += table[i,val-S_i+1] 
            
            end

        end

    end

    return table[size_S, N+1]

end

S = [1,5,10]
N = 12
get_changeways(S, N)



################################################################################
# Problem 2: Rod Cutting
################################################################################

function get_rodvalue(P)
    
    #Length of rod
    n = length(P)

    #Max value of rod lengths 0, 1, 2, ..., n
    V = zeros(n+1) 

    #Loop over lengths
    for ℓ in 1:n
        
        #Initialize a max value
        max_val = -1e10

        #Loop over possible places to cut
        for cut = 1:(ℓ)

            #Check if this is a better place to cut than previous max val
            max_val = max(max_val, P[cut] + V[ℓ - cut + 1])
        
        end

        #Save max val of rod length ℓ
        V[ℓ + 1] = max_val

    end

    V[n+1]

end

P = [1,5,8,9,10,17,17,20]
get_rodvalue(P)



################################################################################
# Problem 3: Knapsack
################################################################################

function knapsack(W::Vector{Int64}, V::Vector{Int64}, C::Int64) 
    #accepts knapsack capacity, weights, and values

    #Number of items
    n = length(W)

    # Columns are capacities 0, 1, ..., C. Rows are set of items up to none, 1, ..., n.
    # i,j entry will be the max val of a knapsack capacity j using items up to i
    # zero capacity and using no items yields zero value. 
    table = zeros(n+1, C+1)

    # Loop over items
    for item in 1:n

        #Loop over capacities 👜
        for 👜 in 1:C

            #Default is value if we don't take this item
            table[item+1,👜+1] = table[item,👜+1]

            #If this item could fit in this bag
            if 👜 >= W[item]

                #Check if we want to include item in bag
                table[item+1,👜+1] = maximum([table[item,👜+1-W[item]] + V[item]; 
                                                table[item,👜+1]])
            
            end

        end

    end

    return table[n+1,C+1]

end


W = [10,20,30]
V = [60,100,120]
C = 50
table = knapsack(W, V, C)


