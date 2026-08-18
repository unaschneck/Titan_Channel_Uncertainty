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

    make_plot = 0;

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
        figure('Name','River Widths')
        hold on;
        
        % Plot river banks
        plot(lon_left, lat_left, 'b.', 'LineWidth', 1.5);
        plot(lon_right, lat_right, 'r.', 'LineWidth', 1.5);
        
        % Plot minimum-width lines
        for i = 1:numel(lat_left)
            j = nearest_right(i);
        
            plot([lon_left(i), lon_right(j)], ...
                 [lat_left(i), lat_right(j)], ...
                 'k-', 'LineWidth', 0.5);

             % Midpoint of the width line
            mid_lon = (lon_left(i) + lon_right(j)) / 2;
            mid_lat = (lat_left(i) + lat_right(j)) / 2;
        
            % Display width value
            if mod(i,10) == 0
            text(mid_lon, mid_lat, sprintf('%.1f', river_width(i)), 'FontSize', 8, 'Color', 'k', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
            end

        end

        % Add 1000 m scale bar
        scale_length = 1000;
        lat_start = lat_left(1);
        lon_start = lon_left(1);
        % Offset from beginning of river
        lat_offset = -0.002;   
        lon_offset =  0.002;  
        x_bar = lon_start + lon_offset;
        y_bar = lat_start + lat_offset;
        % Calculate longitude corresponding to 1000 m at this latitude
        [lat_bar_end, lon_bar_end] = reckon(y_bar, x_bar, scale_length, 90, ellipsoid);
 
        plot([x_bar lon_bar_end],  [y_bar lat_bar_end], 'k-', 'LineWidth', 3);
        text((x_bar + lon_bar_end)/2, y_bar + 0.01 * (max(lat_left) - min(lat_left)), ...
             '1000 m', ...
             'FontSize', 10, ...
             'FontWeight', 'bold', ...
             'Color', 'k', ...
             'HorizontalAlignment', 'center');
        xlabel('Longitude');
        ylabel('Latitude');
        legend('Left bank', 'Right bank', 'Minimum width');
        title('River Width Measurements');
        axis equal;
        grid on;
        hold off;

    end
end