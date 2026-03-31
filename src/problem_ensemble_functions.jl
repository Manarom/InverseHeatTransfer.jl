
"""
    multistart_prob_fun(prob , i , repeat)

Function for solving problems in parallel (can be used for stochastic solver )
```julia
    using Optimization # main interface
    using OptimizationMetaheuristics # PSO optimizer 

    (start, lb, ub)  =  IHT.fill_starting_vectors(parallel_probls) # parallel_probls is precreated inverse problem object
	optp = OptimizationProblem(IHT.discrepancy!, start, parallel_probls, lb = lb, ub = ub , maxiters=pso_iters )
    ensemble_prob = Optimization.EnsembleProblem(
                    optp, 
                    prob_func = multistart_prob_fun
        end
    )

    res_ensemble = solve(ensemble_prob, OptimizationMetaheuristics.PSO(N = 30), EnsembleThreads(), trajectories = 4)

```
"""
function multistart_prob_fun(prob , i , repeat)
    remake(prob , p = deepcopy(prob.p))
end

function regularization_scan(alphas , after_call )
    function prob_func(prob, i, repeat)
        p_new = deepcopy(prob.p)
        α = alphas[i] 
        set_regularization_multiplier!(p_new , α)
        @show α
        return after_call(prob, p = p_new)
    end
    return prob_func
end


function find_l_corner_index(x, y)
    # first and last point 
    p1 = [x[1], y[1]]
    p2 = [x[end], y[end]]
    
    distances = Float64[]
    
    for i in eachindex(x)
        p0 = [x[i], y[i]]
        # distance
        d = abs((p2[2]-p1[2])*p0[1] - (p2[1]-p1[1])*p0[2] + p2[1]*p1[2] - p2[2]*p1[1]) / 
            sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
        push!(distances, d)
    end
    return argmax(distances)
end