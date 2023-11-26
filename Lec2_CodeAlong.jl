############## 0 . machine code ############
# use the @code_native macro to see the machine language that Julia is creating

@code_native 2+2

############# 1. local scope vs. global scope #################

x = 2 # x is created in global scope

# functions, for loops, if statements will create a local scope

function create_z()
    z = 3 # z is created in local scope
end

create_z()
z

for i in 1:10
    w = i + 1
end
w # notice this! w is created in local scope

# Tip 1: Write things inside local scope of function can increase speed since functions get compiled

# example: 
# write in global scope
m = 3
count = 0
@time for i in 1:10^7
    count += m*i
    count -= m*i
end
# write in local scope

function update_count(m)
    count = 0
    for i in 1:10^7
        count += m*i
        count -= m*i
    end
end

@time update_count(3)

# so, always do your calculations inside functions!

############# 2. Types and Structs ############

# typeof() tells you the type of a variable in Julia

# intergers with 64 bits
typeof(1)

# floats with 64 bits
typeof(1.0)

# characters
typeof('a')

# strings
typeof("a")

# vectors
typeof(ones(10))

# matrices
typeof(zeros(2,2))

# 4-dim array of rand
typeof(rand(2,2,2,2))

# you can define type of the variable using :: (this can be used to catch errors if you accidentally not using a right type)
x = 2::Int64
typeof(x)
x = 2::Float64

    # create a new type called MyType
    struct MyType
        # list the fields of that struct
        a
        b
    end

    ex_MyType = MyType(25.0, "number")

    # to access fields, you can just use ".fieldname"
    ex_MyType.a
    ex_MyType.b
    typeof(ex_MyType.a)

    ex2_MyType = MyType("abc", 10)

# Tip 2: Must declare types within your struct in order to increase speed

# example
struct Container_Any
    a # do not declare types
end

struct Container_Float
    a::Float64 # declare types
end

c_any = Container_Any(1.0)
c_float = Container_Float(1.0)

function f(x)
    tmp = x.a
    for i in 1:10^7
        tmp = i + x.a
    end
    tmp
end

@time f(c_any)
@time f(c_float)

