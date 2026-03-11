function [field] = pBCs(field, region, varargin)
%%%%% Applies specified pressure boundary conditions

[ny, nx] = size(field);

switch region
    case "wall"
        field(1,:) = field(2,:);
        field(ny,:) = field(ny-1,:);
    case "inlet"
        field(2:ny-1,1) = field(2:ny-1,2);
    case "outlet"
        field(2:ny-1,nx) = varargin{1};
end