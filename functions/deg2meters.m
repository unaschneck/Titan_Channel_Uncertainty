function [east,north,up] = deg2meters(lat,lon,ellipsoid,lat0,lon0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert distance from degrees to meters (ENU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   lat                 : latitude of points (deg)
%   lon                 : longitude of points (deg)
%   ellipsoid           : ellipsoid to use (e.g., ellipsoid = [2575000,0] for Titan)
%   lat0                : reference latitude (deg)
%   lon0                : reference longitude (deg)
% OUTPUT
%   east                : tangent to Titan's surface pointing towards increasing longitude (m) 
%   north               : tangent to Titan's surface pointing towards increasing latitude (m)
%   up                  : pointing away from surface radially (will decrease or increase due to Titan's curvature relative to the tangent plane) (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    planet_R = ellipsoid.SemimajorAxis;

    lat = deg2rad(lat);
    lon = deg2rad(lon);

    lat0 = deg2rad(lat0);
    lon0 = deg2rad(lon0);

    % spherical to cartesian coordinates
    X = planet_R * cos(lat) .* cos(lon);
    Y = planet_R * cos(lat) .* sin(lon);
    Z = planet_R * sin(lat);

    % Reference point
    X0 = planet_R * cos(lat0) * cos(lon0);
    Y0 = planet_R * cos(lat0) * sin(lon0);
    Z0 = planet_R * sin(lat0);

    % local east, north, up (ENU) basis vectors at starting point
    % ENU coordinates are formed from a plane tangent to the Earth's surface fixed to a specific location
    % https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    east_vec = [-sin(lon0); cos(lon0); 0];
    north_vec = [-sin(lat0)*cos(lon0); -sin(lat0)*sin(lon0); cos(lat0)];
    up_vec = [cos(lat0)*cos(lon0); cos(lat0)*sin(lon0); sin(lat0)];

    % apply basis vectors to points relative to starting position to get
    % coordinates in meters from starting location
    dX = X - X0;
    dY = Y - Y0;
    dZ = Z - Z0;

    east  = east_vec(1)*dX + east_vec(2)*dY + east_vec(3)*dZ;
    north = north_vec(1)*dX + north_vec(2)*dY + north_vec(3)*dZ;
    up    = up_vec(1)*dX + up_vec(2)*dY + up_vec(3)*dZ;

end