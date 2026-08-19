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
%   swath_resolution = meter per pixel at the center of the feature


persistent LAT LON swath seg % make these persistent so counter doesn't have to be reloaded each time which takes forever

if isempty(LAT)

    S = load('counters.mat','lat','lon','swath','seg');

    % swath is a 3D array where 
    %   first dimension corresponds to latitude
    %   second dimension corresponds to longitude
    %   third dimension corresponds to number of observation of that point (up to 10)
    %   the value at each index (i,j,k) corresponds to the SWATH NUMBER that saw the latitude
    %   i and longitude j for observation number k
    % seg is the corresponding array to swath that is also 3D where
    %   first dimension corresponds to latitude
    %   second dimension corresponds to longitude
    %   third dimension corresponds to number of observation of that point (up to 10)
    %   the value at each index (i,j,k) corresponds to the SEGMENT number of the swath (from i,j,k in swath)
    %   that saw the latitude i and longitude j for observation number 

    % Build lat/lon grid
    [LON,LAT] = meshgrid(S.lon(1,:),S.lat(:,1));
  
    swath     = S.swath;      % 5760 x 11520 x 10
    seg       = S.seg;        % 5760 x 11520 x 10

end

make_plot = 1;

% points of the river
feature_lat = feature_lat(:);
feature_lon = feature_lon(:);

% find resolution at center of river
feature_lat_center       = mean(feature_lat(:));
feature_lon_center       = mean(feature_lon(:));
d2feature                = distance(feature_lat_center, feature_lon_center, LAT, LON);  
[~, mindx]               = min(d2feature(:));
[row_center, col_center] = ind2sub(size(LAT), mindx);
lat_POI                  = LAT(row_center,col_center);
lon_POI                  = LON(row_center,col_center);

if col_center < size(LAT, 2) % so doesn't fall of grid
    neighbor_col = col_center + 1;
else
    neighbor_col = col_center - 1;
end
lat_neighbor      = LAT(row_center, neighbor_col);
lon_neighbor      = LON(row_center, neighbor_col);
one_pixel_in_m    = distance(lat_POI,lon_POI,lat_neighbor,lon_neighbor,ellipsoid); 
swath_resolution  = one_pixel_in_m;

disp('...searching for nearest swaths...')
% get associated swath and segment number
swath_all = zeros(numel(feature_lat),size(swath,3),'like',swath);
seg_all   = zeros(numel(feature_lat),size(swath,3),'like',swath);
% find pixels containing river
for k = 1:numel(feature_lat)
    lon_diff = mod(LON - feature_lon(k) + 180, 360) - 180; % wraps around longitude near prime meridian
    dist = lon_diff.^2 + (LAT - feature_lat(k)).^2;
    [~, ind]       = min(dist(:));
    [row,col]      = ind2sub(size(LON),ind);
    swath_all(k,:) = reshape(swath(row,col,:), 1, []);
    seg_all(k,:)   = reshape(seg(row,col,:), 1, []);
end


candidate_pairs = unique([swath_all(:), seg_all(:)], 'rows');

% remove any swaths that don't contain the entire river
swath_number = [];
final_seg = [];
for p = 1:size(candidate_pairs,1)

    this_swath = candidate_pairs(p,1);
    this_seg   = candidate_pairs(p,2);

   if this_seg ~= 1
       continue
   end
    total_covered = false(numel(feature_lat),1);
    
    for pt = 1:numel(feature_lat)
        total_covered(pt) = any(swath_all(pt,:) == this_swath & seg_all(pt,:) == this_seg); 
    end

     if all(total_covered) 
        swath_number = [swath_number; this_swath];
        final_seg = [final_seg; this_seg];
        fprintf('...swath %i (segment %i) includes ALL of the river (%d of %d points)\n',this_swath,this_seg,nnz(total_covered),numel(feature_lat))
     else
        fprintf('...swath %d (segment %i) REMOVED: only %.1f%% of the river (%d of %d points)\n', this_swath, this_seg, 100*nnz(total_covered)/numel(feature_lat), nnz(total_covered), numel(feature_lat))
     end

end




if isempty(swath_number)
    error('No main swath contains the entire river')
end

if ~all(final_seg == 1)
    disp('swath    segment')
    disp([swath_number(:) final_seg(:)])
    error('non-main segment retained')
end

fprintf('...closest swaths: %s\n', strjoin(string(swath_number),', '));

if make_plot
    
    plot_data = swath;
    mask_data = false(size(plot_data));
    for n = 1:numel(swath_number)
        mask_data = mask_data | (plot_data == swath_number(n));
    end

    mask_2D = any(mask_data, 3);

    my_colors = [
    0.90 0.10 0.10   % red
    0.10 0.40 0.90   % blue
    0.10 0.70 0.25   % green
    0.95 0.60 0.05   % orange
    0.60 0.20 0.80   % purple
    0.90 0.20 0.60   % pink
    0.10 0.75 0.75   % cyan
    0.60 0.40 0.10   % brown
    ];
    

    RGB = zeros([size(mask_2D), 3]);
    colored = false(size(mask_2D));

    for n = 1:numel(swath_number)

        m = any(plot_data == swath_number(n), 3);
        m = m & ~colored;

        RGB(:,:,1) = RGB(:,:,1) + m * my_colors(n,1);
        RGB(:,:,2) = RGB(:,:,2) + m * my_colors(n,2);
        RGB(:,:,3) = RGB(:,:,3) + m * my_colors(n,3);

        colored = colored | m;
    end

    figure

    imagesc(LON(1,:), LAT(:,1), RGB, 'AlphaData', 0.8 * mask_2D);

    set(gca,'YDir','normal')
    xlabel('Longitude')
    ylabel('Latitude')
    axis equal tight

    hold on

    % Plot river feature
    plot(feature_lon, feature_lat, '.k', 'MarkerSize', 8)

    % Legend
    for n = 1:numel(swath_number)
        plot(nan,nan,'s', ...
            'MarkerFaceColor', my_colors(n,:), ...
            'MarkerEdgeColor', my_colors(n,:), ...
            'MarkerSize', 8);
    end

    legend(string(swath_number), 'Location','best')


end