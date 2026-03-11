function [] = Results(domain, results)
%%%%% Produce Result Plots

% H2 mass fraction
figure('Name','H_2 Mass Fraction');
contourf(domain.x, domain.y, results.Y_H2, 20, 'LineColor','none'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('H_2 Mass Fraction');
set(gca,'YDir','normal');

% H2O mass fraction
figure('Name','H_2 Mass Fraction');
contourf(domain.x, domain.y, results.Y_H2O, 20, 'LineColor','none'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('H_2O Mass Fraction');
set(gca,'YDir','normal');

% u distribution
figure('Name','u Distribution');
contourf(domain.x, domain.y, results.u, 20, 'LineColor','none'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('u Distribution (m/s)');
set(gca,'YDir','normal');

% v distribution
figure('Name','v Distribution');
contourf(domain.x, domain.y, results.v, 20, 'LineColor','none'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('v Distribution (m/s)');
set(gca,'YDir','normal');

% density
figure('Name','Density');
contourf(domain.x, domain.y, results.rho, 20, 'LineColor','none'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Density (kg/m^3)');
set(gca,'YDir','normal');

% quiver
figure('Name','Flow Direction');
quiver(domain.x, domain.y, results.u, results.v, 'Color', [0.1 0.1 0.7]); % Dark blue arrows
xlabel('x (m)');
ylabel('y (m)');
title('Flow Direction');
axis equal; % Keep aspect ratio to prevent distorted arrows
set(gca,'YDir','normal'); % Ensure y-axis is oriented correctly

% % computing and plotting parabolic distribution
% r = linspace(0, domain.y_L, domain.ny/2)';
% u_r = u(1:domain.ny/2, length(u));
% u_r_analytical = u_inf*1.5*(1 - (r/domain.y_L - 1).^2);
% figure('Name','v Distribution');
% plot(r, u_r, r, u_r_analytical);
% title('Radial u Distribution (m/s)');
% xlabel('r (m)')
% ylabel('u (m/s)')
end