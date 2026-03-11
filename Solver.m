function [results] = Solver(domain, solve, fluid, elec)
%%%%% Solves a pressure-velocity field using projection method
%%% Inputs:
%%% - domain
%%% - solve
%%% - fluid
%%% - elec
%%% Outputs:
%%% - results - all relevant fields


%% Unpacking Inputs
% Solver params
flow_iters = solve.flow_iters;
p_iters_max = solve.p_iters_max;
p_tol = solve.p_tol;
omega_p = solve.omega_p;
omega_i = solve.omega_i;
starve_tol = solve.starve_tol;
dt_safety_factor = solve.dt_safety_factor;

% Geometric Domain
dx = domain.dx;
dx2 = domain.dx2;
dy = domain.dy;
dy2 = domain.dy2;

% Fields
u = domain.u;
v = domain.v;
u_star = domain.u_star;
v_star = domain.v_star;
p = domain.p;
Y_H2 = domain.Y_H2;
Y_H2O = domain.Y_H2O;
i_surface = domain.i_surface;

% BCs
u_inf = fluid.BC.u_inf;
p_out = fluid.BC.p_out;
Y_H2_in = fluid.BC.Y_H2_in;
Y_H2O_in = fluid.BC.Y_H2O_in;

% Indices
i = domain.index.i;
j = domain.index.j;
iplus1 = domain.index.iplus1;
iminus1 = domain.index.iminus1;
jplus1 = domain.index.jplus1;
jminus1 = domain.index.jminus1;

% Fluid params
mu = fluid.mu;
D_H2 = fluid.D_H2;
D_H2O = fluid.D_H2O;
R = fluid.R;

% Electrochemistry params
F = elec.F;
T = elec.T;
M_H2 = elec.M_H2;
M_H2O = elec.M_H2O;
E_rev = elec.E_rev;
E_cell = elec.E_cell;
alpha = elec.alpha;
Resistance = elec.Resistance;
i0 = elec.i0;








%% Applying BCs
u = UBCs(u, "inlet", u_inf);
u = UBCs(u, "wall");
p = pBCs(p, "outlet", p_out);
Y_H2 = speciesBCs(Y_H2, "inlet", Y_H2_in);
Y_H2O = speciesBCs(Y_H2O, "inlet", Y_H2O_in);

rho = p.*LocalMolarMass(Y_H2, Y_H2O, M_H2, M_H2O)/(R*T);
rho_inf = rho(2,1);


%% Calculating dt
C_stab = 0.15; %safety parameter. typically using Explicit Euler, C=0.3-0.5
dt_cfd = C_stab*dx/u_inf;
dt_diffuse = 1*rho_inf/(2*mu)*(1/dx^2 + 1/dy^2)^-1;
dt = min(dt_cfd, dt_diffuse) * dt_safety_factor;


%% Placeholder Arrays
Y_H2_new = Y_H2;
Y_H2O_new = Y_H2O;
mass_total_prev = 0;



%% Solution Loop
for iter = 1:flow_iters

    %% Setup
    % Identifying positive and negative u,v for upwind derivatives
        % Upwind needs +ve for stability 
    u_pos = u(j,i) >= 0;
    u_neg = u(j,i) < 0;
    v_pos = v(j,i) >= 0;
    v_neg = v(j,i) < 0;


    %% Electrochemistry    
    % Partial pressures
    % P_H2_wall = X_H2_wall .* p(1,:); 
    P_H2_wall = (rho(2,:) .* Y_H2(2,:) .* R .* T) ./ M_H2;
    P_H2_ref = 101325;

    % Local reversible potential from Nernst eq.
    E_loc = E_rev + (R*T/(2*F)) * log(max(1e-5, P_H2_wall / P_H2_ref));

    % Butler-Volmer to solve for surface current
    eta = E_loc - E_cell - i_surface * Resistance;
    i_new = i0 * (P_H2_wall / P_H2_ref) .* exp((alpha * F / (R * T)) * eta);
    
    % TEST
    % i_new = 0 * ones(size(i_surface));

    % Current update using under-relaxation for stability
    i_surface = (1 - omega_i) * i_surface + omega_i * i_new;

    
    %% Bottom Wall Flux
    % Computing mass flux out of bottom wall
    mass_flux_H2 = i_surface*M_H2/(2*F);
    starvation_mask = Y_H2(2,:) < starve_tol;
    mass_flux_H2(starvation_mask) = 0;

    % TEST
    % mass_flux_H2(:) = 0;

    v_wall = -mass_flux_H2 ./ rho(2,:);
    v_wall(1) = 0; % for stability

    % Soving for ghost cell Y_H2
    % Pe_wall = v_wall*dy/D_H2; % Peclet number: ratio of advective to diffusive transport
    Y_H2(1, :) = Y_H2(2, :) - (mass_flux_H2 * dy) ./ (rho(2,:) * D_H2);
    Y_H2O(1, :) = 1 - Y_H2(1, :);

    


    %% Species Transport
    % Hydrogen
    dYH2_dx = u_pos .* (Y_H2(j,i) - Y_H2(j,iminus1)) / dx +...
              u_neg .* (Y_H2(j,iplus1) - Y_H2(j,i)) / dx;
    dYH2_dy = v_pos .* (Y_H2(j,i) - Y_H2(jminus1,i)) / dy +...
              v_neg .* (Y_H2(jplus1,i) - Y_H2(j,i)) / dy;

    d2YH2_dx2 = (Y_H2(j,iplus1) - 2*Y_H2(j,i) + Y_H2(j,iminus1)) / dx2;
    d2YH2_dy2 = (Y_H2(jplus1,i) - 2*Y_H2(j,i) + Y_H2(jminus1,i)) / dy2;
    
    Y_H2_new(j,i) = Y_H2(j,i) + dt*(-u(j,i).*dYH2_dx - v(j,i).*dYH2_dy + D_H2*(d2YH2_dx2 + d2YH2_dy2));

    % Water
    dYH2O_dx = u_pos .* (Y_H2O(j,i) - Y_H2O(j,iminus1)) / dx +...
              u_neg .* (Y_H2O(j,iplus1) - Y_H2O(j,i)) / dx;
    dYH2O_dy = v_pos .* (Y_H2O(j,i) - Y_H2O(jminus1,i)) / dy +...
              v_neg .* (Y_H2O(jplus1,i) - Y_H2O(j,i)) / dy;

    d2YH2O_dx2 = (Y_H2O(j,iplus1) - 2*Y_H2O(j,i) + Y_H2O(j,iminus1)) / dx2;
    d2YH2O_dy2 = (Y_H2O(jplus1,i) - 2*Y_H2O(j,i) + Y_H2O(jminus1,i)) / dy2;
    
    Y_H2O_new(j,i) = Y_H2O(j,i) + dt*(-u(j,i).*dYH2O_dx - v(j,i).*dYH2O_dy + D_H2O*(d2YH2O_dx2 + d2YH2O_dy2));

    % Enfoce sum(mass fractions) = 1
    Y_H2_new = max(Y_H2_new, starve_tol);
    Y_H2O_new = max(Y_H2O_new, starve_tol);
    sumY = Y_H2_new + Y_H2O_new;
    Y_H2_new = Y_H2_new ./ sumY;
    Y_H2O_new = Y_H2O_new ./ sumY;

    % Density
    M_bar = LocalMolarMass(Y_H2_new, Y_H2O_new, M_H2, M_H2O);
    rho_new = p.*M_bar/(R*T);


    %% Momentum
    % Convection
    du_dx = u_pos .* (u(j,i) - u(j,iminus1))/dx + ...
            u_neg .* (u(j,iplus1) - u(j,i))/dx;
    du_dy = v_pos .* (u(j,i) - u(jminus1,i))/dy + ...
            v_neg .* (u(jplus1,i) - u(j,i))/dy;
    u_conv = u(j,i).*du_dx + v(j,i).*du_dy;

    dv_dx = u_pos .* (v(j,i) - v(j,iminus1))/dx + ...
            u_neg .* (v(j,iplus1) - v(j,i))/dx;
    dv_dy = v_pos .* (v(j,i) - v(jminus1,i))/dy + ...
            v_neg .* (v(jplus1,i) - v(j,i))/dy;
    v_conv = u(j,i).*dv_dx + v(j,i).*dv_dy;

    % Viscous
    d2udx2 = (u(j,iplus1) - 2*u(j,i) + u(j,iminus1)) / dx^2;
    d2udy2 = (u(jplus1,i) - 2*u(j,i) + u(jminus1,i)) / dy^2;
    d2vdx2 = (v(j,iplus1) - 2*v(j,i) + v(j,iminus1)) / dx^2;
    d2vdy2 = (v(jplus1,i) - 2*v(j,i) + v(jminus1,i)) / dy^2;
    d2udxdy = (u(jplus1,iplus1) - u(jminus1,iplus1) - u(jplus1,iminus1) + u(jminus1,iminus1)) / (4*dx*dy);
    d2vdxdy = (v(jplus1,iplus1) - v(jminus1,iplus1) - v(jplus1,iminus1) + v(jminus1,iminus1)) / (4*dx*dy);

    u_visc = mu./rho_new(j,i) .* (4/3*d2udx2 + d2udy2 + 1/3*d2vdxdy);
    v_visc = mu./rho_new(j,i) .* (4/3*d2vdx2 + d2vdy2 + 1/3*d2udxdy);

    u_star(j,i) = u(j,i) + dt*(-u_conv + u_visc);
    v_star(j,i) = v(j,i) + dt*(-v_conv + v_visc);
     
    % Enforcing BCs on u_star and v_star
    u_star = UBCs(u_star, "inlet", u_inf);
    u_star = UBCs(u_star, "wall");
    u_star = UBCs(u_star, "outlet");
    v_star = UBCs(v_star, "inlet", 0);
    v_star = UBCs(v_star, "topwall");
    v_star = UBCs(v_star, "fluxwall", v_wall);
    v_star = UBCs(v_star, "outlet");

    % Computing intermediate velocity field properties
    dustar_dx = (u_star(j,iplus1) - u_star(j,iminus1)) / (2*dx);
    dvstar_dy = (v_star(jplus1,i) - v_star(jminus1,i)) / (2*dy);
    div_term = ((rho_new(j,i) - rho(j,i))/dt + rho_new(j,i).*(dustar_dx + dvstar_dy)) / dt;


    %% Poisson
    p_diff = p_tol + 1;
    p_iter = 0;
    while p_diff > p_tol && p_iter < p_iters_max
       p_iter = p_iter + 1;
       p_old = p;
        
        p(j,i) = (1-omega_p)*p(j,i) + omega_p/(2/dx2 + 2/dy2).* ...
            ((p(j,iplus1) + p(j,iminus1)) / dx2 + ...
            (p(jplus1,i) + p(jminus1,i)) / dy2 ...
            - div_term);

        % Enforcing p BC's
        p = pBCs(p, "inlet");
        p = pBCs(p, "outlet", p_out);
        p = pBCs(p, "wall");

        % Computing difference 
        p_diff = max(max(abs(p - p_old)));
    end

    % Updating states
    dpdx = (p(j,iplus1) - p(j,iminus1)) / (2*dx);
    dpdy = (p(jplus1,i) - p(jminus1,i)) / (2*dy);
    u(j,i) = u_star(j,i) - (dt ./ rho_new(j,i)) .* dpdx;
    v(j,i) = v_star(j,i) - (dt ./ rho_new(j,i)) .* dpdy;
    rho = rho_new;
    Y_H2 = Y_H2_new;
    Y_H2O = Y_H2O_new;

    % Enforce BCs after correction
    u = UBCs(u, "inlet", u_inf);
    u = UBCs(u, "wall");
    u = UBCs(u, "outlet");
    v = UBCs(v, "inlet", 0);
    v = UBCs(v, "topwall");
    v = UBCs(v, "fluxwall", v_wall);
    v = UBCs(v, "outlet");
    Y_H2 = speciesBCs(Y_H2, "topwall");
    Y_H2 = speciesBCs(Y_H2, "inlet", Y_H2_in);
    Y_H2 = speciesBCs(Y_H2, "outlet");
    Y_H2O = speciesBCs(Y_H2O, "topwall");
    Y_H2O = speciesBCs(Y_H2O, "inlet", Y_H2O_in);
    Y_H2O = speciesBCs(Y_H2O, "outlet");
    rho = speciesBCs(rho, "outlet");
    rho = speciesBCs(rho, "inlet", rho_inf);
    rho = speciesBCs(rho, "wall");


    %% Diagnostics
    % Iteration number
    disp("Iter No. : " + iter + " / " + flow_iters + " | (" + iter/flow_iters*100 + "%)")
    
    % Pressure convergence
    disp("Pressure iters: " + p_iter)

    % Wall velocity
    disp("v_wall: " + v_wall(1) + " " + v_wall(2) + " " +v_wall(3))

    % % Midchannel u/u_inf
    % disp("u_midchannel/u_inf: " + u(ceil(domain.ny/2),ceil(domain.nx/2))/u_inf)

    % Conversavtion of Mass
    mass_total = sum(rho(:)); % compare before/after step
    mass_change = mass_total - mass_total_prev;
    fprintf('Mass change this step = %.3e\n', mass_change);
    mass_total_prev = mass_total;

    % Face Flux
    flux_left  = sum((rho(:,1) .* u(:,1)) * dy);      
    flux_right = sum((rho(:,end) .* u(:,end)) * dy);  
    flux_bottom = sum((rho(1,:) .* v(1,:)) * dx);     
    flux_top    = sum((rho(end,:) .* v(end,:)) * dx); 
    disp(flux_right + " " + flux_left + " " + flux_top + " " + flux_bottom)
end

%% Packing Results
results.u = u;
results.v = v;
results.p = p;
results.rho = rho;
results.Y_H2 = Y_H2;
results.Y_H2O = Y_H2O;
results.i = i_surface;
results.avg_i = mean(i_surface);
