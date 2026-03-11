function [field] = speciesBCs(field, region, varargin)
%%%%% Applies specified species boundary conditions

[ny, nx] = size(field);

switch region
    case "bottomwall"
        field(1,:) = 0.01; % for testing
    case "topwall"
        field(ny,:) = field(ny-1,:);
    case "wall"
        field(ny,:) = field(ny-1,:);
        field(1,:) = field(2,:);
    case "inlet"
        field(2:ny-1,1) = varargin{1};
    case "outlet"
        field(2:ny-1,nx) = field(2:ny-1,nx-1);
end