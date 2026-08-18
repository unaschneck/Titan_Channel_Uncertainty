function my_sin = calc_sinuosity(x,y,ellipsoid)
% function calculates the sinusoity of the channel where Sin > 1 always
% and larger Sin describes more windy path
% INPUTS:
%       x         = longitude in degrees
%       y         = latitude in degrees
%       ellipsoid = structure containing planet's semimajor axis in meters
% OUTPUTS:
%       my_sin    = sinuosity (alongpath distance / straight path distance)

        R = ellipsoid.SemimajorAxis; % m

        % Convert coordinates to radians
        lon = deg2rad(x);
        lat = deg2rad(y);

        % cumulative distance along the river path
        dLon = diff(lon);
        dLat = diff(lat);
        c =  haversine(lat(1:end-1), lat(2:end), dLat, dLon);
        d_river_path = R * sum(c);

         % distance between endpoints as the crow flies (straight line distance)
        dlon_total = lon(end) - lon(1);
        dlat_total = lat(end) - lat(1);

        c_straight = haversine(lat(1), lat(end), dlat_total, dlon_total);

        d_straight_path = R * c_straight;

        % Sinuosity
        my_sin = d_river_path / d_straight_path;

    function c = haversine(lat_start, lat_end, dLat, dLon)

        a = sin(dLat/2).^2 + cos(lat_start) .* cos(lat_end) .* sin(dLon/2).^2;
        % Great circle distance
        c = 2 .* atan2(sqrt(a), sqrt(1-a));

    end


end