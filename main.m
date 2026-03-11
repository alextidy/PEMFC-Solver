clc
clear
close all

%%%%% Coupled PEMFC Channel Solver
%%% Solver Architecture:
%%% - Fix cell voltage, solve for current distribution
%%% - E_cell is determined using Butler-Volmer
%%% - Distribution of i at the cathode 

%%% Notes: 
%%% - j=2 is the bottom wall, j=nyr is the top
%%% - One ghost cell exists at each boundary

%%% Stabel pipeflow conditions: 
%%% - x_L = 0.5 
%%% - y_L = 0.1
%%% - nxr = 20 
%%% - nyr = 20
%%% - u_inf = 0.05 m/s
%%% - p_out = 101325 Pa

%%% IMPORTANT: Does not model GDL dynamics

%% Constants
% Fluid Constants
fluid.mu = 1.8e-5;
fluid.D_H2 = 1e-6;
fluid.D_H2O = 1e-4;
fluid.R = 8.314;

% Electrochemical Constants
elec.F = 96485;
elec.T = 273 + 25;
elec.M_H2 = 0.002;
elec.M_H2O = 0.018;
elec.E_rev = 0.0;
elec.alpha = 0.5; % Transfer coeff
elec.Resistance = 1e-3; % Total fuel cell resistance
elec.i0 = 1;




%% Options 
% Numerical Options
solve.flow_iters = 10e3;
solve.p_iters_max = 10000;
solve.p_tol = 1e-8;
solve.omega_p = 1.0; % Above 1 is overrelaxation, below is underrelaxation
solve.omega_i = 0.1;
solve.starve_tol = 1e-5;
solve.dt_safety_factor = 1e-1; % Recommended 1e-1 to 1e-3

% Fluid Options
fluid.BC.u_inf = 0.01;
fluid.BC.p_out = 101325;
fluid.BC.Y_H2_in = 0.8;
fluid.BC.Y_H2O_in = 1 - fluid.BC.Y_H2_in;

% Domain Options
domainset.x_L = 0.50;
domainset.y_L = 0.10;
domainset.nxr = 20; % cells in x direction, excluding ghost
domainset.nyr = 20; % cells in y direcition, excluding ghost

% Volatge Options
i_iters = 30;
E_cell_array = linspace(-0.0, -0.8, i_iters);
i_array = zeros(i_iters, domainset.nxr+2);
avg_i_array = zeros(i_iters, 1);



%% Running Simulation
for i = 1:i_iters

    % Diagnostics


    % Reinitialise Domain
    domain = DomainGen(domainset, fluid.BC);

    % New voltage
    elec.E_cell = E_cell_array(i);

    % Solve Domain with specified BCs
    results = Solver(domain, solve, fluid, elec);

    % Store currents
    i_array(i,:) = results.i;
    avg_i_array(i) = results.avg_i;

end


%% Results Processing
% Single run
Results(domain, results);

% Multirun
BuildPolar(E_cell_array, avg_i_array)