
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
        #if isa(p_new , SingleInverseProblem)
        #    p_new.α[] = α
        #elseif isa(p_new , ParallelInverseProblems)
            for p in p_new.problems
                p.α[] = α
            end
        #end
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
## ENSEMBLE AVERAGING
#=
        #optimizer = OptimizationOptimJL.ParticleSwarm()
        #Optimization.EnsembleProblem
        #ensemble_prob = Optimization.EnsembleProblem(optp, safety_copy = true)
        #res = solve(optp, optimizer , use_initial = true)      
        
        ensemble_prob = EnsembleProblem(
                        optp, 
                        prob_func = (prob, i, repeat) -> begin
                # Создаем глубокую копию задачи для каждого потока
                # Это гарантирует, что in-place изменения не перемешаются
                remake(prob, p = deepcopy(prob.p))
            end
            )
        res_ensemble = solve(ensemble_prob, PSO(N = 30), EnsembleThreads(), trajectories = 4)

        best_sol = res_ensemble[argmin([s.objective for s in res_ensemble])]

        lam_coeffs_number = length(IHT.coeffs(inv_probl.optimizable.λ))
        res.u[lam_coeffs_number + 1 : end]/2700






## AI generated L-curve methods  need to be checked!!
error("Unchecked")

 using Optimization, OptimizationMetaheuristics



ensemble_prob = EnsembleProblem(optp, prob_func = prob_func)

# Запускаем параллельно (каждое ядро считает свой вариант альфа)
res_ensemble = solve(ensemble_prob, PSO(N = 40), EnsembleThreads(), trajectories = length(alphas))


# Собираем логарифмические координаты
# x_coords: log(||residual||^2)
# y_coords: log(||regularization||^2)
x_coords = [log10(sol.p.covariance(sol.p.residual, sol.p.t_grid.dt)) for sol in res_ensemble]
y_coords = [log10(sum(finite_difference_regularization_loss(ov) for ov in sol.p.optimizable)) for sol in res_ensemble]


function find_l_corner_index(x, y)
    # Координаты начала и конца кривой
    p1 = [x[1], y[1]]
    p2 = [x[end], y[end]]
    
    distances = Float64[]
    for i in 1:length(x)
        p0 = [x[i], y[i]]
        # Расстояние от точки p0 до прямой p1-p2
        d = abs((p2[2]-p1[2])*p0[1] - (p2[1]-p1[1])*p0[2] + p2[1]*p1[2] - p2[2]*p1[1]) / 
            sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
        push!(distances, d)
    end
    return argmax(distances)
end

best_idx = find_l_corner_index(x_coords, y_coords)
best_alpha = alphas[best_idx]
println("Оптимальное значение альфа: ", best_alpha)

=#