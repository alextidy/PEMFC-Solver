function [field] = UBCs(field, region, varargin)
%%%%% Applies specified velocity boundary conditions

[ny, nx] = size(field);

switch region
    case "wall"
        field(1,:) = 0;
        field(ny,:) = 0;
    case "topwall"
        field(ny,:) = 0;
    case "fluxwall"
        field(1,:) = varargin{1};
    case "inlet"
        field(2:ny-1,1) = varargin{1};
    case "outlet"
        field(2:ny-1,nx) = field(2:ny-1,nx-1);
end