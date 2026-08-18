function [lat_interpolated,lon_interpolated] = even_geodesic_spacing(lat,lon,num_points,ellipsoid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate line to have num_pts points that are evenly spaced in geodesic distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   lat                 : latitude of points on original line
%   lon                 : longitude of points on original line
%   num_points          : numeber of desired points on interpolated line
%   ellipsoid           : ellipsoid to use (e.g., ellipsoid = [2575000,0] for Titan)
% OUTPUT
%   lat_interpolated    : latitude of new points interpolated to even spacing 
%   lon_interpolated    : longitude of new points interpolated to even spacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computes the geodesic distance between every consecutive pair
lon1 = lon(1:end-1);
lon2 = lon(2:end);

lat1 = lat(1:end-1);
lat2 = lat(2:end);

% Segment distances on Titan ellipsoid
len_path = distance(lat1, lon1, lat2, lon2, ellipsoid);

% fraction of path length
new_pts = linspace(0,1,num_points); 

% Normalize lengths
arclen = len_path/sum(len_path);

% cumulative arc length
cumarc = [0;cumsum(arclen)];

% determine which segment each interpolation point belongs to
[~,tbins] = histc(new_pts,cumarc); 
  
% fix boundary conditions (aka first point belongs to first segment and last point to last segment)
tbins((tbins <= 0) | (new_pts <= 0)) = 1;
tbins((tbins >= numel(lat)) | (new_pts >= 1)) = numel(lat) - 1;
  
% compute local interpolation fraction
s = (new_pts' - cumarc(tbins))./arclen(tbins);

% compute segment distance and azimuth
[my_length,my_azimuth] = distance(lat(tbins),lon(tbins),lat(tbins+1),lon(tbins+1),ellipsoid);

% compute new interpolated points
[lat_interpolated,lon_interpolated] = reckon(lat(tbins),lon(tbins),my_length.*s,my_azimuth,ellipsoid);


end