
clearvars
%% initial settings
DEAFULT_TMEASURMENT_NAME = "Tmeasured.txt";
DEAFULT_HEATING_REGIME_NAME = "heating_regime.txt";
% initial data
thermocouple_locations = [3,6.8];% measurements thermocouples locations
Tmeasured_inds = [1,3,4]; % indices in Tmeasured data first  - for time, other for thermocouple locations
thickness = 7.0;
timestep = 1;
control_thermocouple_coordinate = 0.2;
coordinate_step  = 0.05;
min_time = 120;
max_time = 370;
density = 2720;
% initial quess and constraints
lam = [19, 17, 13, 10, 9, 8, 7,7]';
 % lb = 0:numel(lam(:)) - 1;
 % ub = 22 - 2*lb; % lambda upper constrants
 % lb = 15.5 - 1.7*lb; % lambda lower constrants
ub = 22*lam./lam;
lb = 5*lam./lam;
% optimization settings
alfa = 1e-4; % Tikhonov regularization factor
SWARM_ITERATIONS_NUMBER = 1;% number of PSO reruns
IS_SWARM_PLOT = true; % show optimization process 
OPTIMIZER = 2; % optimizer, 1 - NelderMead, 2 - particle swarm
% regularization function
%reg_fun = @(X)transpose(X(:))*X(:)/length(X); %  alfa*I
reg_fun = @(X)transpose(diff(X(:)))*diff(X(:))/length(X);% finite difference

%% loading input data
%cur_file_folder = fullfile(cur_file_folder,"testing_data","lam_inversion");
cur_file_folder = uigetdir();
if isnumeric(cur_file_folder)
    cur_file_folder = cd;%fullfile(cur_file_folder,"testing_data","flux_inversion2");
else
    cd(cur_file_folder)
end
assert(isfolder(cur_file_folder), "Unable to load data from " + cur_file_folder)
% regime_file = fullfile(cur_file_folder, "test_heating_regime.txt"); % heating regime loading
% measure_file = fullfile(cur_file_folder,"test_Tmeasured.txt"); % file with measured temperatures
regime_file = fullfile(cur_file_folder, "heating_regime.txt"); % heating regime loading
measure_file = fullfile(cur_file_folder,"Tmeasured.txt"); % file with measured temperatures

try
    Tmeasured = OneDHeatTransfer.make_unique(OneDHeatTransfer.clear_nans(readmatrix(measure_file)));
    heating_regime =OneDHeatTransfer.make_unique(OneDHeatTransfer.clear_nans(readmatrix(regime_file)));
catch ex
    disp(ex.message)
    return
end
Tmeasured = trim_timedata(Tmeasured,min_time,max_time); %
heating_regime = trim_timedata(heating_regime,min_time,max_time);
%q_dwn = OneDHeatTransfer.eval_radiative_flux(Tmeasured(:,end),20,0.2); % dwn heat flux calculation
q_dwn = @(T)-OneDHeatTransfer.eval_radiative_flux(T,20,0.2);
Tmeasured = Tmeasured(:,Tmeasured_inds);
%thermocouple_locations = [0.7,3,6.8];



%% specifying thermal conductivity

flag = true(size(lam));

cp_starting_vector = [ -0.000000000217660, 0.000001155312137, -0.002192746146536, 1.85569440318715, 650.053870289652000];
cp_optimizable_parameters_flag = [false,false,false,false,false];
%% optimization problem formulation

for iii = 1:SWARM_ITERATIONS_NUMBER
    problems(1,iii) = OneDHeatTransfer(Tmeasured,heating_regime,density, ...
        "thickness",thickness,"dt",timestep, ...
        "control_thermocouple_coordinate",control_thermocouple_coordinate,"q_dwn", ...
         q_dwn,"thermocouple_coordinates",thermocouple_locations,"dh",coordinate_step, ...
        "default_basis","BernsteinBasis");
    % filling starting vector
        problems(iii).set_lam_parametric_function(lam,flag)
        problems(iii).set_cp_parametric_function(cp_starting_vector,cp_optimizable_parameters_flag)
end

%% plotting methods
ax1 = problems(1).plot_data(data="measured");
problems(1).plot_data(ax1,data="heating regime",hold="on");
legend(ax1,["Measured" "Heating regime"])
%% plotting quantities
% ax_lam = get_next_ax();
% ax_dldT = get_next_ax();
% 
% problems(1).plot_quantity(ax_lam,quantity = "lam",hold="on")
% problems(1).plot_quantity(ax_dldT,quantity = "dldT",hold="on");

%% solving direct problem
ax_distr = get_next_ax();
for p =problems
    solve(p)
end
p.plot_data(ax_distr,data ="fitted_temperatures",hold="on");
%% solving the optimization problem

optimization_starting_vector = lam;

fval = zeros(SWARM_ITERATIONS_NUMBER,1);
swarm_res = zeros(length(lb),SWARM_ITERATIONS_NUMBER);

switch OPTIMIZER
    case 1
        if IS_SWARM_PLOT && SWARM_ITERATIONS_NUMBER
            options = optimset('TolX',1e-4,'TolF',1e-4,'PlotFcns',@optimplotfval);
        else
            options = optimset('TolX',1e-4,'TolF',1e-4);
        end
    case 2
        if IS_SWARM_PLOT && SWARM_ITERATIONS_NUMBER
            options = optimoptions("particleswarm","Display","iter","PlotFcn","pswplotbestf");
        else
            options = optimoptions("particleswarm","Display","iter");
        end
end
if SWARM_ITERATIONS_NUMBER > 1
    parfor iii = 1:SWARM_ITERATIONS_NUMBER
        p = problems(iii);
        fun_i = @(X) p.discrepancy(X) +  sum(X(:) < lb(:)) + sum(X(:) > ub(:)) + alfa*(reg_fun(X)); % sum(problem.lam(problem.Theating(:)) < 0) +
        
        switch OPTIMIZER
            case 1
                [xval,fval_cur,exitflag,output] = fminsearch(fun_i,optimization_starting_vector,options);
            case 2
                [xval,fval_cur,exitflag,output] = particleswarm(fun_i,numel(optimization_starting_vector),lb,ub,options);
        end
        fval(iii) = fval_cur;
        swarm_res(:,iii) = xval;
    end
else
    for iii = 1:SWARM_ITERATIONS_NUMBER
        p = problems(iii);
        fun_i = @(X) p.discrepancy(X) +  sum(X(:) < lb(:)) + sum(X(:) > ub(:)) + alfa*(reg_fun(X)); % sum(problem.lam(problem.Theating(:)) < 0) +
        
        switch OPTIMIZER
            case 1
                [xval,fval_cur,exitflag,output] = fminsearch(fun_i,optimization_starting_vector,options);
            case 2
                [xval,fval_cur,exitflag,output] = particleswarm(fun_i,numel(optimization_starting_vector),lb,ub,options);
        end
        fval(iii) = fval_cur;
        swarm_res(:,iii) = xval;
    end
end
%
%% searching for best fit 
[~,min_ind] = min(fval);
optimization_starting_vector = swarm_res(:,min_ind);
%% refit the problem 
problem = problems(min_ind);
problem.refresh_parameters(optimization_starting_vector)
fun = @(X) problem.discrepancy(X) +  sum(X(:) < lb(:)) + sum(X(:) > ub(:)) + alfa*(reg_fun(X));
%optimization_starting_vector = problem.fill_starting_vector();
options2 = optimset('TolX',1e-6,'TolF',1e-6,'PlotFcns',@optimplotfval);
[xval_local,fval_local,exitflag,output] = fminsearch(fun,optimization_starting_vector,options2);

%% plotting solution
% problem.plot_data("data","residual");
% [~,lam_fitted] = problem.plot_quantity(quantity = "lam");
% [~,Cp_fitted] = problem.plot_quantity(quantity = "Cp");
% [~,alpha_fitted] = problem.plot_quantity(quantity = "alpha");
% problem.plot_data("data","distribution_over_coordinate");
% problem.plot_data("data","dwn_heat_flux");
% problem.plot_data("data","upper_heat_flux");
% problem.plot_data("data","supplied_heat_flux");
% problem.plot_data("data","fitted_temperatures");
%% 
problem.plot_data("data","residual");
%% comparing to data 

lam_real_fun = OneDHeatTransfer.material_properties_function("RBSN","lam");
ax_f = problem.plot_quantity(quantity = "lam");
data = problem.get_physical_property("quantity","lam");
T_data = data(:,1);
lam_inv = data(:,2);
lam_real = lam_real_fun(T_data);

lower_bound = problem.lam_polynomial.poly_eval_unnorm(lb,T_data);
upper_bound = problem.lam_polynomial.poly_eval_unnorm(ub,T_data);

ln = ax_f.Children(1);
ln.LineWidth = 3;
hold(ax_f,"on")
plot(ax_f, T_data,lam_real,"r", LineWidth=3)
plot(ax_f, T_data,lower_bound,"g", LineWidth=3)
plot(ax_f, T_data,upper_bound,"m", LineWidth=3)

hold(ax_f,"off")
legend(ax_f,["inverse", "real","lower","upper"])
save("problem.mat","problem")