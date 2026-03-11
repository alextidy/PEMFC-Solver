function [domain] = DomainGen(domainset, BC)
%%%%% Sets up a problem domain with ghost cells to solve
%%% Assumes a linearly spaced rectangular grid

% Size constraints
domain.x_L = domainset.x_L;
domain.y_L = domainset.y_L;
domain.nxr = domainset.nxr;
domain.nyr = domainset.nyr;
domain.nx = domainset.nxr + 2;
domain.ny = domainset.nyr + 2;
domain.dx = domain.x_L / domain.nxr;
domain.dx2 = domain.dx^2;
domain.dy = domain.y_L / domain.nyr;
domain.dy2 = domain.dy^2;

% Indexing arrays
domain.index.i = 2:domain.nxr+1;
domain.index.j = 2:domain.nyr+1;
domain.index.iplus1 = 3:domain.nxr+2;
domain.index.iminus1 = 1:domain.nxr;
domain.index.jplus1 = 3:domain.nyr+2;
domain.index.jminus1 = 1:domain.nyr;

% Fields
domain.u = zeros(domain.ny,domain.nx);
domain.v = zeros(domain.ny,domain.nx);
domain.u_star = zeros(domain.ny,domain.nx);
domain.v_star = zeros(domain.ny,domain.nx);
domain.p = ones(domain.ny,domain.nx) * BC.p_out;
domain.Y_H2 = ones(domain.ny,domain.nx) * BC.Y_H2_in;
domain.Y_H2O = ones(domain.ny,domain.nx) * BC.Y_H2O_in;
domain.i_surface = ones(1,domain.nx) * 100;
[domain.x, domain.y] = meshgrid(linspace(0,domain.x_L,domain.nx), linspace(0,domain.y_L,domain.ny));
end