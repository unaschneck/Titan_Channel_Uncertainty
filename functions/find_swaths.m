function [swath_number,swath_resolution] = find_swaths(feature_lat,feature_lon,ellipsoid)
% identify the swaths that contain input coordinates
%
% INPUT:
%   feature_lat = latitude coordinates of river feature
%   feature_lon = longitude coordinates of river feature
%   ellipsoid   = ellipsoid representation of Titan for distance calculations
%
% OUTPUT:
%   swath_number = swath IDs containing at least one feature point
%   swath_resolution = meter prt one pixel at the center of the feature


persistent LAT LON swath_name sz % make these persistent so counter doesn't have to be reloaded each time which takes forever

if isempty(LAT)

    S = load('counters.mat','lat','lon','swath','seg');

    % Build lat/lon grid
    [LON,LAT] = meshgrid(S.lon(1,:),S.lat(:,1));
  
    swath_name = S.swath;      % 5760 x 11520 x 10
    sz         = size(S.lat);

end

make_plot = 0;

feature_lat = feature_lat(:);
feature_lon = feature_lon(:);


% find resolution at center of river
feature_lat_center = mean(feature_lat(:));
feature_lon_center = mean(feature_lon(:));
d2feature = distance(feature_lat_center, feature_lon_center, LAT, LON);  
[~, mindx] = min(d2feature(:));
[row_center, col_center] = ind2sub(size(LAT), mindx);
lat_POI = LAT(row_center,col_center);
lon_POI = LON(row_center,col_center);
if col_center < size(LAT, 2) % so doesn't fall of grid
    neighbor_col = col_center + 1;
else
    neighbor_col = col_center - 1;
end
lat_neighbor = LAT(row_center, neighbor_col);
lon_neighbor = LON(row_center, neighbor_col);
one_pixel_in_m = distance(lat_POI,lon_POI,lat_neighbor,lon_neighbor,ellipsoid); 
swath_resolution = one_pixel_in_m;

disp('....searching for nearest swaths...')

% knnsearch finds the nearest grid cell, grid spacing is regular
min_idx = min((LON(:)-feature_lon).^2 + (LAT(:)-feature_lat).^2); %knnsearch([LAT(:) LON(:)], [feature_lat feature_lon]);
[row, col] = ind2sub(sz, min_idx);


% Extract swath IDs at each nearest grid cell
num_layers = size(swath_name, 3);
% convert back to linear index
lin_idx = sub2ind(sz, row, col);

swaths_at_points = zeros(numel(lin_idx), num_layers, class(swath_name));

% collect all swaths that overlap with the river
for k = 1:num_layers
    layer = swath_name(:, :, k);
    swaths_at_points(:, k) = layer(lin_idx);
end
swath_number = unique(swaths_at_points(:));
swath_number = swath_number(swath_number ~= 0);

if make_plot
    data = swath_name(:,:,1);
    mask = ismember(data, swath_number);
    RGB = zeros([size(data), 3]);
    colors = [
        0.90 0.10 0.10   % red
        0.10 0.40 0.90   % blue
        0.10 0.70 0.25   % green
        0.95 0.60 0.05   % orange
        0.60 0.20 0.80   % purple
    ];
    for n = 1:numel(swath_number)
        m = data == swath_number(n);
    
        RGB(:,:,1) = RGB(:,:,1) + m * colors(n,1);
        RGB(:,:,2) = RGB(:,:,2) + m * colors(n,2);
        RGB(:,:,3) = RGB(:,:,3) + m * colors(n,3);
    end
    figure
    h = image(LON(1,:), LAT(:,1), RGB);
    set(gca,'YDir','normal')
    
    % Make everything not selected transparent
    h.AlphaData = mask;
    
    xlabel('Longitude')
    ylabel('Latitude')
    % Legend
    hold on
    plot(feature_lon,feature_lat,'.k')
    
    for n = 1:numel(swath_number)
        plot(nan,nan,'s', ...
            'MarkerFaceColor',colors(n,:), ...
            'MarkerEdgeColor',colors(n,:), ...
            'MarkerSize',8);
    end
    
    legend(string(swath_number), 'Location','best')
end

% Keep only swaths that contain all of the feature
swath_number_filtered = zeros(1, numel(swath_number), class(swath_number));
num_found = 0;
num_feature_points = size(swaths_at_points, 1);

for i = 1:numel(swath_number)

    % feature points in this swath
    inside_track = any(swaths_at_points == swath_number(i), 2);
    % fraction of feature contained within swath
    fraction_feature = nnz(inside_track) / num_feature_points;

    if fraction_feature == 1.0

        num_found = num_found + 1;
        swath_number_filtered(num_found) = swath_number(i);

        fprintf('...swath %d includes ALL of the river (%d of %d points)\n',swath_number(i), nnz(inside_track),  num_feature_points);

    else

        fprintf('...swath %d REMOVED: only %.1f%% of the river (%d of %d points)\n', swath_number(i), 100 * fraction_feature, nnz(inside_track), num_feature_points);

    end

end

swath_number = swath_number_filtered(1:num_found);

fprintf('...closest swaths: %s\n', strjoin(string(swath_number),', '));

end