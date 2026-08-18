function [node_spacing,az_spacing] = find_spacing(x,y,ellipsoid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find spacing between consecutive points (e.g., on the centerline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   x            : x coordinates of line
%   y            : y coordinates of line
%   ellipsoid    : ellipsoid of planet for calculating distance
% OUTPUT
%   node_spacing : spacing between points on a line (in dimensions of x,y)
%   az_spacing   : azimuthal spacing between points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    pt1 = [y(1:numel(x)-1), x(1:numel(x)-1)];
    pt2 = [y(2:numel(x)), x(2:numel(x))];
    [node_spacing,az_spacing] = distance(pt1,pt2,ellipsoid);
   

end