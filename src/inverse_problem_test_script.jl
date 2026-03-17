	
    #using Revise
    include(joinpath(@__DIR__,"InverseHeatTransfer.jl"))
    using Plots, BenchmarkTools, Test
    using CSV, Tables , Optimization,OptimizationOptimJL
    using OptimizationNLopt
    using OptimizationMetaheuristics
    import .InverseHeatTransfer as IHT
    # working_folder = raw"E:\JULIA\JULIA_DEPOT\dev\InverseHeatTransfer.jl\test\test_data\property_inversion_ansys_new\25ks"
    working_folder = raw"D:\JuliaDepoth\dev\InverseHeatTransfer.jl\test\test_data\property_inversion_ansys_new\25ks"
    hr_file = filter(f->  contains(f,"Tmeasured") , readdir(working_folder))
	isempty(hr_file) && error("folder must constain file with  Tmeasured")
	T_measured = CSV.read(joinpath(working_folder,hr_file[]), Tables.matrix)

    tmax = 150.0
    measured_inds = (2, 3, 4 , 5 , 6, 11)
    thermocouples_locations = [1e-3*(0.3 + (i - 2)) for i in measured_inds]   # 1e-3 * [0.3 , 1.3 , 4.3 , 9.3]
    time_data = T_measured[:,1]
	temperatures = T_measured[: , collect(measured_inds)]
    fl = time_data .<= tmax
    time_data = time_data[fl]
    temperatures = temperatures[fl, :]
    plot(time_data, temperatures)
    lam_T_range = [extrema(temperatures)...]


	initial_distribution = lam_T_range[1]
    thickness = 9.3*1e-3
    
    optimize_cp_flag = true
    basis_degree = 5
    xpoints_number = 200
    tpoints_number = 2000
    Cp_default_vector = (2700*900.0 , 2700*1000.0 , 2700*1000.0 , 2700*1000.0 , 2700*1000.0 , 2700*1000.0)
    Cp_default_vector = ntuple(_-> 2700*1000.0 , basis_degree)

    C_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(Cp_default_vector) , 
        xmin = lam_T_range[1], xmax = lam_T_range[end])
    C = IHT.OptimizableVariable(C_poly, 
        lb = 2700*300.0, ub = 2700*1500.0, 
        flag = optimize_cp_flag)

	λ_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_-> 1.0 , basis_degree)), xmin = lam_T_range[1], xmax = lam_T_range[end])
	λ = IHT.OptimizableVariable(λ_poly, flag = true , lb = 3.0, ub=25.0)
	dλdT_poly =  IHT.PolynomialWrappers.derivative(λ_poly)
	dλdT = IHT.OptimizableVariable(dλdT_poly)

    props_T_range = range(lam_T_range[1],lam_T_range[2],100)
    plot(props_T_range, C.(props_T_range))
    plot(props_T_range, λ.(props_T_range))

    inv_probl = IHT.SingleInverseProblem(time_data, 
        temperatures,
        initial_distribution,
        thermocouples_locations,
        C,λ, dλdT, thickness,
        xpoints_number, tpoints_number; 
        regularization = IHT.NoRegularization(),
        covariance = IHT.AR1Covariance(2.071 , 1.0))


        #@benchmark IHT.regularization_loss($inv_probl)
        #@benchmark IHT.constraints_loss($inv_probl)
        #@benchmark IHT.modify!($(inv_probl.optimizable.λ) , $[5.0, 4.0, 8.0] )
    
        IHT.discrepancy!( [ 21.0, 15.0, 10.2 , 2700*1e3 , 2700*1e3 , 2700*1e3] , inv_probl)
        @benchmark IHT.discrepancy!( $[ 21.0, 15.0, 10.2 , 2700*1e3 , 2700*1e3 , 2700*1e3] , $inv_probl)

        plot(inv_probl.residual)

        (x0, lb, ub) = IHT.fill_starting_vectors(inv_probl)
	    problems = (inv_probl , deepcopy(inv_probl) , deepcopy(inv_probl))
        
        parprob = IHT.ParallelInverseProblems(problems...)

        @benchmark IHT.discrepancy!($x0, $parprob)
        @benchmark IHT.discrepancy!($x0, $inv_probl)

        optp = OptimizationProblem(IHT.discrepancy!, x0, parprob, lb = lb, ub = ub  ,  maxiters=100)
	    #optimizer = PSO(N = 10*2*basis_degree,  options = Metaheuristics.Options(parallel_evaluation = true))
        optimizer = OptimizationNLopt.NLopt.LN_COBYLA()


        #Optimization.EnsembleProblem
        #ensemble_prob = Optimization.EnsembleProblem(optp, safety_copy = true)
       #res = solve(optp, optimizer )         #
      #=  
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
=#