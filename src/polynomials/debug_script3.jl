using StaticArrays, Polynomials, Plots
include("PolynomialWrappers.jl")

#V = PolynomialWrappers.VanderMatrix(SVector((1.0,2.3,4.5,6.7,8.9)),PolynomialWrappers.StandPolyWrapper{10,Float64})
a = (1.0,1.0,1.0,2.3,4.5,6.7)
p = PolynomialWrappers.TrigPolyWrapper(a)
p_der = PolynomialWrappers.derivative(p)
x = -1.0:1e-2:1.0
ppp = Polynomials.fit(Polynomial,x,p.(x),20)

plot(x,p.(x))
plot!(x,ppp.(x))

ppp_der = Polynomials.derivative(ppp)
plot(x,ppp_der.(x))
plot!(x,p_der.(x))


function ai_eval_scaled_poly(ch::ChebPoly{N,T}, x::S, a, b) where {N,T,S}
    R = promote_type(T, S)
    
    # Обработка пустых или константных полиномов
    N == 0 && return zero(R)
    N == 1 && return R(ch.coeffs[1])
    
    # Масштабирование в [-1, 1]
    ξ = R(2 * (x - a) / (b - a) - 1)
    two_ξ = 2 * ξ
    
    # Алгоритм Кленшоу (Clenshaw)
    # Инициализируем b_{n+1} и b_{n+2} как 0
    b_next2 = zero(R) # b_{k+2}
    b_next1 = zero(R) # b_{k+1}
    
    @inbounds for k in N:-1:2
        # b_k = c_k + 2ξ * b_{k+1} - b_{k+2}
        b_curr = R(ch.coeffs[k]) + two_ξ * b_next1 - b_next2
        b_next2 = b_next1
        b_next1 = b_curr
    end
    
    # Финальный шаг для T_0: c_1 + ξ * b_2 - b_3
    return R(ch.coeffs[1]) + ξ * b_next1 - b_next2
end


function ai_eval_scaled_poly(leg::LegPoly{N,T}, x::S, a, b) where {N,T,S}
    R = promote_type(T, S)
    
    # Защита для низких степеней
    N == 0 && return zero(R)
    N == 1 && return R(leg.coeffs[1])
    
    # Масштабирование в [-1, 1]
    ξ = R(2 * (x - a) / (b - a) - 1)
    
    b_next2 = zero(R) # b_{k+2}
    b_next1 = zero(R) # b_{k+1}
    
    # Идем от N до 2 (обратный ход Кленшоу)
    @inbounds for k in N:-1:2
        n = k - 1 # Порядок полинома P_n
        
        # Коэффициенты для b_k = c_k + α_k*ξ*b_{k+1} - β_{k+1}*b_{k+2}
        # Для Лежандра: α_n = (2n+1)/(n+1), β_{n+1} = (n+1)/(n+2)
        α = R(2n + 1) / R(n + 1)
        
        # Важно: β берется для следующего шага рекурсии
        # На шаге k мы вычисляем b_k, используя b_{k+1} и b_{k+2}
        # Множитель при b_{k+2} для Лежандра это (n+2)/(n+1) * β_{n+2}
        # Упрощенная форма Кленшоу для Лежандра:
        β_val = R(n + 1) / R(n + 2) # Это β_{k+1}
        
        # Корректная формула шага:
        b_curr = R(leg.coeffs[k]) + α * ξ * b_next1 - (R(k) / R(k+1)) * b_next2
        
        b_next2 = b_next1
        b_next1 = b_curr
    end
    
    # Финальный шаг (P_0 = 1): c_1 + α_0*ξ*b_2 - β_1*b_3
    # Для n=0: α_0 = 1, β_1 = 1/2
    return R(leg.coeffs[1]) + ξ * b_next1 - R(0.5) * b_next2
end