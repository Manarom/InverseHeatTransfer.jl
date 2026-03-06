	
    #using Revise
    include(joinpath(@__DIR__,"InverseHeatTransfer.jl"))
    using Plots, BenchmarkTools, Test
    using CSV, Tables
    import .InverseHeatTransfer as IHT
    working_folder = raw"E:\JULIA\JULIA_DEPOT\dev\InverseHeatTransfer.jl\test\test_data\property_inversion_ansys_new\25ks"
    hr_file = filter(f->  contains(f,"Tmeasured") , readdir(working_folder))
	isempty(hr_file) && error("folder must constain file with  Tmeasured")
	T_measured = CSV.read(joinpath(working_folder,hr_file[]), Tables.matrix)

    measured_inds = (2,3,6,11)
    thermocouples_locations = [0.3,1.3,4.3,9.3]
    time_data = T_measured[:,1]
	temperatures = T_measured[: , collect(measured_inds)]
	initial_distribution = 20.0
    thickness = 9.3

    lam_T_range = [20,1500]
    basis_degree = 3
    xpoints_number = 200
    tpoints_number = 2000

    C = IHT.OptimizableVariable(IHT.ScaledPolynomial(IHT.BernsteinSymPoly((2700*900.0 , 2700*1000.0 , 2700*1000.0)) , xmin = lam_T_range[1], xmax = lam_T_range[end]), flag = true)
	λ_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->1.0 , basis_degree)), xmin = lam_T_range[1], xmax = lam_T_range[end])
	λ = IHT.OptimizableVariable(λ_poly, flag = true , lb = 3.0,ub=25.0)
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
        xpoints_number, tpoints_number)


        is_optim = Base.Fix2(isa, IHT.OptimizableVariable)
        #=optimizable = (; is_optim(t[2]) for t in zip((:λ,:C),(λ, C))]
        for z in zip((:λ,:C),(λ, C))
            @show typeof(z) 
            @show is_optim(z[2])

        end=#

        nt = (; (k => v for (k, v) in zip((:λ,:C),(λ, C)) if is_optim(v))...)

        for ov in nt
            @show typeof(ov)
        end


        IHT.update_all_optimizables!(inv_probl, [2.0,2.0,2.0, 2700*1300, 2700*1300, 2700*1350])

        plot(props_T_range, inv_probl.optimizable.λ.(props_T_range))
        plot(props_T_range, inv_probl.optimizable.dλdT.(props_T_range))

        IHT.coeffs(λ)

        IHT.fill_starting_vector(inv_probl)