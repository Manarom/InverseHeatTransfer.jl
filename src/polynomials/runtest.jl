using StaticArrays, Test, BenchmarkTools, Plots,LinearAlgebra, Test
include("PolynomialWrappers.jl")
import Polynomials, LegendrePolynomials, SpecialPolynomials
import .PolynomialWrappers as PW
show_title(st) = begin 
    println("-------------------------------------")
    println(st)
    println("-------------------------------------")
end
show_subtitle(st,n = 1) = begin
    sh = mapreduce(_-> "==",*,1 : n)
    println("$sh > $(string(st))")
end

@testset "Polynomials" begin 
    #begin 
        a = (0.1,-0.5,0.4,-1.0)
        N = length(a)
        show_title("Testing polynomials evaluation:")
        @test Polynomials.Polynomial(a)(0.5) ≈ PW.StandPoly(a)(0.5)
        @test  Polynomials.ChebyshevT(a)(0.5) ≈ PW.ChebPoly(a)(0.5)
        @test  LegendrePolynomials.Pl(0.5 , 3) ≈ PW.LegPoly((0.0, 0.0, 0.0, 1.0))(0.5)
        @test  SpecialPolynomials.Bernstein{N - 1}(a)(0.5) ≈ PW.BernsteinPoly(a)(0.5)

        a = [1.0 , 22.0 , 3.0 , 4.56 , 3.51]
        N = length(a)
        x = range(-1.0,1.0,10)
        show_subtitle("Testing polynomial derivative coefficients and benchmarking the evaluation speed")
        for (TestPolyType,CheckPolyType) in zip((SpecialPolynomials.Bernstein{N - 1},
                                SpecialPolynomials.Chebyshev,
                                SpecialPolynomials.Legendre ), (PW.BernsteinSymPoly, 
                            PW.ChebPoly, PW.LegPoly))
                show_subtitle("Testing  $(CheckPolyType)", 2)
                b_test = TestPolyType(a)
                b_der_test = Polynomials.derivative(b_test)
                b_pw = CheckPolyType(a)
                b_der = PW.derivative(b_pw)
                show_subtitle("$(TestPolyType)(x) evaluation: ", 3)
                @btime $b_test(0.5)
                show_subtitle(" $(CheckPolyType)(x) evaluation: ", 3)
                @btime $b_pw(0.5)
                
                show_subtitle("derivative( $(TestPolyType) ) evaluation: ",3)
                @btime $Polynomials.derivative($b_test)
                show_subtitle("derivative( $(CheckPolyType) ) evaluation: ",3)
                @btime $PW.derivative($b_pw)

                if ~(TestPolyType <: SpecialPolynomials.Bernstein)
                    @test norm(PW.coeffs(b_der) .- b_der_test.coeffs)≈ 0 atol = 1e-12
                else # BErnstein basis scaled from 0 to  1
                    @test norm(PW.coeffs(b_der) .- b_der_test.coeffs./2)≈ 0 atol = 1e-12
                end

        end

        x = collect(range(-1.0,1.0,100))
        a = (1.0,22.0,3.0,4.56,3.51)
        show_title("Default scaling polynomial derivative evaluation test:")
        for CheckPolyType in (PW.BernsteinSymPoly, PW.StandPoly, 
                                    PW.ChebPoly, PW.LegPoly, PW.TrigPoly)

            show_subtitle("Cheking the derivative for $CheckPolyType")
            pw_poly = CheckPolyType(a)
            pw_poly_vals = pw_poly.(x)
            is_trig = CheckPolyType <: PW.TrigPoly

            s_fit = Polynomials.fit(Polynomials.Polynomial, x, pw_poly_vals, is_trig ? 20 : PW.poly_degree(pw_poly)) 
    
            s_der = Polynomials.derivative(s_fit)
            pw_poly_der = PW.derivative(pw_poly)
            @test norm(s_der.(x) .- pw_poly_der.(x)) ≈ 0 atol= is_trig ? 1e-10 : 1e-12 # need to chenge the accuracy for trig because of bad fitting
        end


        #Scaled polynomial derivatives
        x_unscaled = collect(200.0:0.1:2000.0)
        a = (1.0,22.0,3.0,4.56,3.51)
        (x_scaled,x_min,x_max) = PW.scale_x_to_ξ(x_unscaled)

        p_scaled = PW.ScaledPolynomial(PW.ChebPoly, a,  x_min, x_max)
        y_values = p_scaled.(x_unscaled)
        p_scaled_der = PW.derivative(p_scaled)
        pp = Polynomials.fit(Polynomials.Polynomial,x_unscaled, y_values, length(a) - 1)
        y_der_test = Polynomials.derivative(pp).(x_unscaled)

        plot(x_unscaled,y_der_test)
        plot!(x_unscaled, p_scaled_der.(x_unscaled))

        show_title("Scaled polynomial derivative evaluation test ")
        for CheckPolyType in (PW.BernsteinSymPoly, PW.StandPoly, PW.ChebPoly , PW.LegPoly)
            show_subtitle("Cheking the derivative for $CheckPolyType")
            p_scaled = PW.ScaledPolynomial(CheckPolyType, a,  x_min, x_max)
            y_values = p_scaled.(x_unscaled)
            @test norm(y_values .- CheckPolyType(a).(x_scaled)) ≈ 0.0 atol =1e-12
            # 
            pp = Polynomials.fit(Polynomials.Polynomial,x_unscaled, y_values, length(a) - 1)
            y_der_test = Polynomials.derivative(pp).(x_unscaled)

            p_scaled_der = PW.derivative(p_scaled)

            @test norm( p_scaled_der.(x_unscaled) .- y_der_test ) ≈ 0.0 atol = 1e-12

        end
        println("-------------------------------------")

end
