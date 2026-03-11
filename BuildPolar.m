function [] = BuildPolar(E_cell_array, avg_i_array)
%%%%% Build anode polar
%%% Small voltage range due to solver stability constraints
%%% - Likely fixable with implicit methods

% Standardising voltage and calculating power
E_cell_array = 1 - abs(E_cell_array);
avg_power = E_cell_array.*avg_i_array';
avg_power = avg_power ./ max(avg_power(~isnan(avg_power)));

% Plotting average current against E_cell
figure('Name', 'Anode Polarization Curve');
hold on
plot(avg_i_array, E_cell_array, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r')
plot(avg_i_array, avg_power, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b')
grid on;
xlabel('Average Current Density (A/m^2)');
ylabel('Driving Force |E_{cell}| (V)');
title('Anode Half-Cell Polarization Curve');
legend("Voltage", "Normalised Power Density", "Location","best")
end