function [M_bar] = LocalMolarMass(Y_H2, Y_H20, M_H2, M_H2O)
%%%%% Compute local molar mass from hyrogen and water

M_bar = 1 ./ (Y_H2/M_H2 + Y_H20/M_H2O);
end