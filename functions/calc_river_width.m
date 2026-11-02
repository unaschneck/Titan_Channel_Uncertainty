function river_width = calc_river_width(lat_left,lon_left,lat_right,lon_right,ellipsoid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the minimum distance between banks to determine the river width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%       lat_left    : left bank latitude (interpolated to even spacing) [1 x n]
%       lon_left    : left bank longitude (interpolated to even spacing) [1 x n]
%       lat_right   : right bank latitude (interpolated to even spacing) [1 x n]
%       lon_right   : left bank longitude (interpolated to even spacing) [1 x n]
%   OUTPUT
%       river_width : width of river moving downstream [1 x n]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    make_plot = 1;

    river_width = NaN(size(lat_left));
    nearest_right = NaN(size(lat_left));


    for i = 1:numel(lat_left)
        min_distance = inf;
    
        for j = 1:numel(lat_right)
    
            pairwise_distance = distance(lat_left(i), lon_left(i), ...
                lat_right(j), lon_right(j), ...
                ellipsoid);
    
            if pairwise_distance < min_distance
                min_distance = pairwise_distance;
                nearest_right(i) = j;
            end
        end
    
        river_width(i) = min_distance;
    end
    
    if make_plot

        figure('Name', 'River Widths');
        
        geobasemap satellite;
        
        % Plot river banks
        geoplot(lat_left, lon_left, 'b.', 'LineWidth', 1.5);
        hold on;
        geoplot(lat_right, lon_right, 'r.', 'LineWidth', 1.5);
        
        % Plot minimum-width lines
        for i = 1:numel(lat_left)
        
            j = nearest_right(i);
        
            geoplot([lat_left(i), lat_right(j)], ...
                    [lon_left(i), lon_right(j)], ...
                    'k-', 'LineWidth', 0.5);
        
            % Midpoint of the width line
            mid_lon = (lon_left(i) + lon_right(j)) / 2;
            mid_lat = (lat_left(i) + lat_right(j)) / 2;
        
            % Display width value every 10 points
            if mod(i, 100) == 0
                text(mid_lat, mid_lon, ...
                    sprintf('%.1f', river_width(i)), ...
                    'FontSize', 8, ...
                    'Color', 'k', ...
                    'BackgroundColor', 'w', ...
                    'HorizontalAlignment', 'center');
            end
        
        end
        

    end
end