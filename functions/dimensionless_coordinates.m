function [x_dimensionless,y_dimensionless] = dimensionless_coordinates(lat,lon,avgwidth,ellipsoid,lat0,lon0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% converts lat/lon coordinates to dimensionless cartesian coordinates measured
% east-west and north-south normalized to the average width of the river
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%      lat             : latitude in geodesic coordinates (deg)
%      lon             : longitude in geodesic coordinates (deg)
%      avgwidth        : average river width, used to make final coordinate system dimensionless
%      ellipsoid       : ellipsoid of planet
%   OUTPUT
%      x_dimensionless : east-west Cartesian coordinates normalized to average river width
%      y_dimensionless : north-south Cartesian coordinates normalized to average river width
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % convert lat,lon in deg to meters relative to first point
   [east_west,north_south,~] = deg2meters(lat,lon,ellipsoid,lat0,lon0);

    % normalize by avg width
    x_dimensionless = east_west./avgwidth;
    y_dimensionless = north_south./avgwidth;


end


