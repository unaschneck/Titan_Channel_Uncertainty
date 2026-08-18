function [cum_dis,az] = cumulative_distance_downstream(lat_center,lon_center,ellipsoid)

    lat_center = lat_center(:);
    lon_center = lon_center(:);

    lat1 = lat_center(1:end-1);
    lon1 = lon_center(1:end-1);
    
    lat2 = lat_center(2:end);
    lon2 = lon_center(2:end);
    
    
    
    [dpath, az] = distance(lat1, lon1, lat2, lon2, ellipsoid);

    cum_dis = [0; cumsum(dpath)];




end