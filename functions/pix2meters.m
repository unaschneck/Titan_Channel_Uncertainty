function [east,north] = pix2meters(col_pix,row_pix,sar_file,ellipsoid,lat0,lon0)

line   = row_pix;                                                          % row index
sample = col_pix;                                                          % column index

[lat,lon] = bidr_proj_inv_v2(sample,line,sar_file);                        
[east,north,~] = deg2meters(lat,lon,ellipsoid,lat0,lon0);

end

