function [row_pix,col_pix] = transform_coordinates(lat,lon,sar_file)
% project feature coordinates into image bidr coordinates using Alex's 
% external function that scrapes metadata from the sbdr files
%
% INPUT:
%   lat = latitude coordinates of river feature
%   lon = longitude coordinates of river feature
%   sar_file = sar file name containing the feature data
%
% OUTPUT:
%   col_pix   = East-West pixel index of feature (column)
%   row_pix   = North-South pixel index of feature (row)


[line,sample] = bidr_proj_v2(lat,lon,sar_file); % function found on Alex's W drive

row_pix = line;     % row index
col_pix = sample;   % column index


end