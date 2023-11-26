############ 1. comment ############
# this is a one-line comment

#=
    this is a multi-line comment
=#


########### 2. Math #############

# Julia is just a fancy calculator
2+2

10^4 + 2*3/100

# suppressing and printing output

# semi-colon can suppres output
2+2;
ans # this will give the last-line output
ans = 3;
3+4
ans # calculation will override the value of ans

# print function will write output to the REPO
print(3*4);

# the @show macro will print the code and the answer to the terminal
    # macro is basically a function of the codeline
@show 5^3

########### 3. variables ############
 
# we can save results of our calculation as variables

x = 2*3
@show x

my_number = 225

β   # \beta
γ   # \gamma
🐟  # \:fish
🚗  # \:car

# there are also special operations for updating variables

x = 1
x += 1
x -= 1
x /= 2
x *= 3

########## 4. arrays, vectors, and matrices #############

# an array is a collection of elements of the same type
# a vector is a 1-dim array
v1 = rand(4) # this create 4 random numbers in range 0 to 1
v2 = collect(1:2:20)
v3 = collect(range(start = 1, stop = 100, length = 199))

# a matrix is a 2-dim array
A = rand(6,4)
B = zeros(4,2)
C = ones(2,3)
D = ones(10)

# matrix operations

    # multiplication
    x = rand(4)
    Z = rand(5,4)
    Z*x
    # transpose
    Z'
    # inverse
    Y = rand(4,4)
    inv(Y)

# broadcasting

A = ones(5,5)
log(A) # this will not work
log.(A)
A.*3

# ctrl+c will attemp to stop your code line from running
# or we can hit the trash can simbol

rand(5,6,7,8) # we can define any dim as we want

############## 5. packages #############

# type the ] in REPO will open the package interface, for example, we add package "Distributions" by typing "] add Distributions"

# load the package into the current environment
using Distributions
using Plots
using Random

############# 6. For loops #############

# construct a dataset x_{t+1} = x_{t} + ϵ where ϵ ~ N(0, 1)
x_vec = zeros(100)

Random.seed!(234) # set up the random seed
for i in 2:100 # loop over 2nd element to 100th element
    x_vec[i] = x_vec[i-1] + rand(Normal())
end

x_vec
plot(x_vec);


y_vec = zeros(100)

Random.seed!(123) # set up a different random seed
for i in 2:100 # loop over 2nd element to 100th element
    y_vec[i] = y_vec[i-1] + rand(Normal())
end

y_vec
plot(y_vec)

################ 7. defining functions #############

function OLS(y, X)
    X'
    2+2
    inv(X'*X)*X'*y # the last line will be the return value
end

X = rand(Normal(), 100)
ϵ = rand(Normal(), 100)
y = 2X + ϵ

OLS(y, X)

############## 8. Booleans and Conditionals ##############
true
false

2 == 3

5 > 1

7 <= 10

8 != 8


x = 2
if x > 5
    print("x is greater than 5")
elseif x > 3
    print("x is greater than 3 but leq to 5")
else
    print("x is leq than 3")
end

# While loops

i = 1

while i < 100
    i += 1
end

i

# run another julia file
include("other_file.jl")

# make my own arrays
x = [2 3 4;
     5 6 7;
     8 9 10]
A = [rand(2,2) rand(2,2)]
A[1,4]
A[2,2] = 2
A[1,:] = [0 0 0 0]
A
