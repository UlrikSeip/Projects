#julia integrator class
#essentially an ode solver with several different integration methods
#object oriented syntax in julia:
"""
type MyType
a::Int64
b::Float64
end

x = MyType(3, 4)

x.a

---------------------
function double(x::MyType)
    x.a *= 2
end
"""


type ODEsolver
    #insert stuff