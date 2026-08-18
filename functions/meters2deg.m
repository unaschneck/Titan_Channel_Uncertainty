function [lat,lon] = meters2deg(east,north,up,ellipsoid,lat0,lon0)

planet_R = ellipsoid.SemimajorAxis;

% Convert reference location to radians
lat0 = deg2rad(lat0);
lon0 = deg2rad(lon0);

% Reference ECEF position
X0 = planet_R * cos(lat0) * cos(lon0);
Y0 = planet_R * cos(lat0) * sin(lon0);
Z0 = planet_R * sin(lat0);

% Local ENU basis vectors at reference point
east_vec = [-sin(lon0); cos(lon0); 0];
north_vec = [-sin(lat0)*cos(lon0); -sin(lat0)*sin(lon0); cos(lat0)];
up_vec = [cos(lat0)*cos(lon0);  cos(lat0)*sin(lon0); sin(lat0)];

% Convert local ENU displacement back to ECEF displacement
dX = east_vec(1)*east +  north_vec(1)*north + up_vec(1)*up;
dY = east_vec(2)*east + north_vec(2)*north + up_vec(2)*up;
dZ = east_vec(3)*east + north_vec(3)*north + up_vec(3)*up;

% ECEF coordinates
X = X0 + dX;
Y = Y0 + dY;
Z = Z0 + dZ;

% ECEF -> spherical latitude/longitude
lon = atan2(Y,X);
lat = atan2(Z,sqrt(X.^2 + Y.^2));

% Convert radians -> degrees
lat = rad2deg(lat);
lon = rad2deg(lon);

end