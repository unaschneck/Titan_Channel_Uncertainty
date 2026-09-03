clc
clear
close all
% restoredefaultpath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT TO RUN THE EARTH COMPARISON OF RIVER FEATURE DEPENDANCE
% ON RESOLUTION
% created 8/2026 by U. Schneck (ugschneck@gmail.com)
% This script will loop through the three Earth river analogues in 
% Northern Alaska, Northern Quebec, and Northeastern Australia to determine
% how feature geometry is affected by changes to resolution.
% Figures produced:
%   (1) Width as a function of downstream distance and resolution
%   (2) Width as a funciton of resolution
%   (3) Channel brightness as a function of resolution
% 
% This script expects mapped channel data located in 
% data/EARTH_ANALOGUES/ which contain a folder for each  river at high 
% resolution (30 m/pixel), medium  (175 m/pixels), and low  resolution
% (350 m/pixel). These maps were created by F. Verdi 2023 by degrading
% the resolution of original SAR images at 30 m/pixel through downsampling
% and adding speckle noise. 
% Ex data structure
%   |--- Alaska
%           |-- Arc                  : ArcGIS project where mapping was performed
%           |
%           |-- geometry             : geometry and radar_data structures (post-analysis), including width and brightness
%           |     |--- high-res
%           |     |--- medium-res
%           |     |--- low-res
%           |
%           |-- rasters              : SAR data (.tif)
%           |     |--- high-res
%           |     |--- medium-res
%           |     |--- low-res
%           |
%           |-- shapefiles            : mapped centerline and banks (.shp)
%                 |--- high-res
%                 |--- medium-res
%                 |--- low-res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_subdivisions = 10;

data_dir = fullfile(pwd,'data','EARTH_ANALOGUES');

results_dir = 'results';
if ~isfolder(results_dir)
    mkdir(results_dir);
end


RESOLUTIONS = [350 175 30]; % m/pixel
EARTH_RIVERS = {'Alaska','Australia','Canada'}; 
RES_FILES = {'high-res','medium-res','low-res'};

% run analysis on each river at each resolution
% each resolution has a different map

ellipsoid = make_ellipsoid('Earth');                                                 % simple sphere with radius of Titan for geographic projection
num_points = 1000;                                                                   % number of points for the final mapped channels (will be evenly spaced)
max_num_subdivisions = 10;                                                           % maximum number of subdivisions to consider for cross and along channel brightness measurement

avg_width = NaN(numel(EARTH_RIVERS),numel(RESOLUTIONS));
std_width = NaN(numel(EARTH_RIVERS),numel(RESOLUTIONS));

res_width_fig = figure('Name','Res vs Width');
ax1 = gca;
xlabel('Resolution (m/pixel)')
ylabel('Width (m)')
axis([0 400 0 800])
grid on;
hold on

mycolor = orderedcolors("dye");


for i = 1:numel(EARTH_RIVERS)

    disp('==================================================================')
    fprintf('River in %s:\n',upper(EARTH_RIVERS{i}));

    for j = 1:numel(RES_FILES)

            shp_data_loc = fullfile(data_dir,EARTH_RIVERS{i},'shapefiles',RES_FILES{j});

            % extract centerline and bank data from shapefile
            [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = extract_from_shp(shp_data_loc,'Earth');
            if isempty(x_leftedge)
                warning('no available shapefile for %s at resolution %i m/pix\n',EARTH_RIVERS{i},RESOLUTIONS(j))
                continue
            end

            % make sure the lines are mapped in the same orientation downstream
            [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = fix_mapping_order(x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge);

            % interpolate the centerline to be evenly spaced
            [lat_center,lon_center] = even_geodesic_spacing(y_center,x_center,num_points,ellipsoid);
            % interpolate the left bank
            [lat_left,lon_left] = even_geodesic_spacing(y_leftedge,x_leftedge,num_points,ellipsoid);
            % interpolate the right bank
            [lat_right,lon_right] = even_geodesic_spacing(y_rightedge,x_rightedge,num_points,ellipsoid);

            % calculate river width along channel
            river_width = calc_river_width(lat_left,lon_left,lat_right,lon_right,ellipsoid);

            avg_width(i,j) = mean(river_width,'omitmissing');                       % mean width of river 
            std_width(i,j) = std(river_width,'omitmissing');                        % std of river width

            fprintf('Average width (resolution = %i) = %0.2f ± %0.2f\n',RESOLUTIONS(j),avg_width(i,j),std_width(i,j))

            % normalize coordinates of centerline to average river width
            [xd_center,yd_center] = dimensionless_coordinates(lat_center,lon_center,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));
            % normalize coordinates of left bank to average river width
            [xd_left,yd_left] = dimensionless_coordinates(lat_left,lon_left,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));
            % normalize coordinates of right bank to average river width
            [xd_right,yd_right] = dimensionless_coordinates(lat_right,lon_right,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));

            % find spacing between points normalized by average river width
            [node_spacing,az_spacing] = find_spacing(lon_center,lat_center,ellipsoid);
            % normalize spacing to the average river width
            node_d_spacing = node_spacing./mean(river_width,'omitmissing');

            % find width of channel at each point the centerline
            width_centerline = width_along_centerline(lon_center,lat_center,lon_left,lat_left,river_width,ellipsoid);

            % find cumulative distance downstream
            downstream_distance = cumulative_distance_downstream(lat_center,lon_center,ellipsoid);
            % convert left bank coordinates into distance from centerline
            downstream_left = calc_dis_edge(lon_center,lat_center,lon_left,lat_left,ellipsoid);

            path_downstream_fig = figure('Name','Channel path and downstream width');

            subplot(2,1,1)
            scatter(lon_center,lat_center,10,1:length(lat_center),"filled")
            colormap("winter")
            grid on;
            xlabel('lon (deg)')
            ylabel('lat (deg)')
            title(sprintf('%s river at %i m/pix',EARTH_RIVERS{i},RESOLUTIONS(j)))

            subplot(2,1,2)
            std_top = river_width + std_width(i,j);
            std_bot = river_width - std_width(i,j);
            xx2 = [downstream_distance', fliplr(downstream_distance')];
            std_shade = [std_top', fliplr(std_bot')];
            fill(xx2,std_shade,[0.7 0.7 0.7],'EdgeColor','none')
            hold on;
            scatter(downstream_distance,river_width,10,1:length(downstream_distance),'filled')
            colormap("winter")
            xlabel('downstream distance (km)')
            ylabel('river width (m)')
            grid on;
            ylim([0 max(river_width)+100])
            drawnow;

            radar_data_loc = fullfile(data_dir,EARTH_RIVERS{i},'rasters',[RES_FILES{j} '.tif']);
           
            % extract brightness data of channel
            [x_feature,y_feature,radar_vals,brightness_vals] = extract_earth_RADAR(radar_data_loc,lon_left,lat_left,lon_right,lat_right);

            % transform into meters relative to starting position in local east, north, up coordinates (ENU)
           [x_center, y_center]           = deg2meters(lat_center,lon_center,ellipsoid,lat_center(1),lon_center(1));
           [x_left, y_left]               = deg2meters(lat_left,lon_left, ellipsoid,lat_center(1),lon_center(1));
           [x_right, y_right]             = deg2meters(lat_right,lon_right,ellipsoid,lat_center(1),lon_center(1));

           % calculate downstream brightness
           
           [downstream_brightness_distance,brightness_profile] = calc_brightness_profile_downstream(x_feature,y_feature,radar_vals,brightness_vals,x_left,y_left,x_right,y_right,downstream_distance,num_subdivisions)
           % 
           %  % calculate cross stream brightness and contrast
           % 



    end

    errorbar(ax1,RESOLUTIONS,avg_width(i,:),std_width(i,:),'-s','Color',mycolor(i,:),'LineWidth',2,'DisplayName',EARTH_RIVERS{i})

end

legend(ax1,'show')

%saveas(channel_fig,fullfile(results_dir,['EARTH_channel_path.png']))
