function width_centerline = width_along_centerline(x_center,y_center,x_left,y_left,width,ellipsoid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the width associated with each point on the centerline (primarily for making
% maps look better)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   x_center            : longitude of points on centerline (deg)
%   y_center            : latitude of points on centerline (deg)
%   x_left              : longitude of points on left bank (deg)
%   y_left              : latitude of points on left bank (deg)
%   x_right             : longitude of points on right bank (deg)
%   y_right             : latitude of points on right bank (deg)
%   width               : array of river widths (m)
%   ellipsoid           : ellipsoid to use (e.g., ellipsoid = [2575000,0] for Titan)
% OUTPUT
%   width_centerline    : array of rivers width associated with each point downstream (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    width_centerline = zeros(numel(x_center),1);

    for i = 1:numel(x_center)
        edge_dis = distance(y_center(i)*ones(numel(x_left),1),x_center(i)*ones(numel(x_left),1),y_left, x_left, ellipsoid);

        [~, edge_ind] = min(edge_dis);
        width_centerline(i) = width(edge_ind);
    end

end