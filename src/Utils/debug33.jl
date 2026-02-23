using FunctionWrappers, BenchmarkTools


struct A{T}
    f::FunctionWrappers.FunctionWrapper{T,Tuple{T}}
    A(f, ::Type{T}) where T= new{T}(f) 
end

struct B{F}
    f::F 
end

function loop_fun(a)
    s = a(0.5)
    for _ in 1:1000
        s += a(0.5)
    end
    return s
end

(p::A{T})(x::T) where T = p.f(x)
(p::B)(x::T) where T= p.f(x)::T

a = A(sin, Float64)
b = B(sin)
a(0.5)
b(0.5)

@benchmark a(0.5)
@benchmark b(0.5)

@benchmark loop_fun($a)
@benchmark loop_fun($b)

@code_warntype loop_fun(a)