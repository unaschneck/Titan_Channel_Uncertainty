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
axis([0 400 0 1400])
grid on;
hold on;


mycolor = orderedcolors("dye");


for i = 1:numel(EARTH_RIVERS)

    disp('==================================================================')
    fprintf('River in %s:\n',upper(EARTH_RIVERS{i}));

    for j = 1:numel(RES_FILES)

            shp_data_loc = fullfile(data_dir,EARTH_RIVERS{i},'Una','shapefiles',RES_FILES{j});

            % extract centerline and bank data from shapefile
            [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = extract_from_shp(shp_data_loc,'Earth');
            if isempty(x_leftedge) | isempty(x_center) | isempty(x_rightedge)
                error('no available shapefile for %s at resolution %i m/pix\n',EARTH_RIVERS{i},RESOLUTIONS(j))
            end
            % 
            % figure;
            % plot(x_center,y_center,'-k')
            % hold on
            % plot(x_leftedge,y_leftedge,'-r')
            % plot(x_rightedge,y_rightedge,'-b')

            % make sure the lines are mapped in the same orientation downstream
            [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = fix_mapping_order(x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge);
    
            % figure;
            % plot(x_center,y_center,'-k')
            % hold on
            % plot(x_leftedge,y_leftedge,'-r')
            % plot(x_rightedge,y_rightedge,'-b')

            % interpolate the centerline to be evenly spaced
            [lat_center,lon_center] = even_geodesic_spacing(y_center,x_center,num_points,ellipsoid);
            % interpolate the left bank
            [lat_left,lon_left] = even_geodesic_spacing(y_leftedge,x_leftedge,num_points,ellipsoid);
            % interpolate the right bank
            [lat_right,lon_right] = even_geodesic_spacing(y_rightedge,x_rightedge,num_points,ellipsoid);

            % figure;
            % plot(lon_center,lat_center,'-k')
            % hold on
            % plot(lon_left,lat_left,'-r')
            % plot(lon_right,lat_right,'-b')
            % 
            % figure;
            % geoplot(lat_center,lon_center,'-r')
            % hold on
            % geoplot(lat_left,lon_left,'-b')
            % geoplot(lat_right,lon_right,'-k')
            % geobasemap satellite

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

            saveas(path_downstream_fig,fullfile(results_dir,['EARTH_channelpath_' EARTH_RIVERS{i} '_' num2str(RESOLUTIONS(j)) '.png']))


            all_width{i,j} = river_width;



    end

    errorbar(ax1,RESOLUTIONS,avg_width(i,:),std_width(i,:),'-s','Color',mycolor(i,:),'LineWidth',2,'DisplayName',EARTH_RIVERS{i})

end

legend(ax1,'show')

%saveas(res_width_fig,fullfile(results_dir,'EARTH_channel_width_v_res.png'))


mycolor1 = hex2rgb('#ab3428');
mycolor2 = hex2rgb('#2d728f');
close all

figure('Name','Width Error')
hold on

for i = 1:3

    true_width = all_width{i,1}(:);
    high_width = all_width{i,2}(:);
    low_width  = all_width{i,3}(:);

    % Valid paired observations
    valid_high = ~isnan(true_width) & ~isnan(high_width);
    valid_low  = ~isnan(true_width) & ~isnan(low_width);

    diff_high = true_width(valid_high) - high_width(valid_high);
    diff_low  = true_width(valid_low) - low_width(valid_low);

    histogram(diff_high,'BinEdges',-1005:50:1005,'FaceColor',mycolor1,'DisplayName','175 m/pixel')
    histogram(diff_low,'BinEdges',-1005:50:1005,'FaceColor',mycolor2,'DisplayName','350 m/pixel (Cassini)')
    xline(mean(diff_high(:),'omitmissing'),'--','Color',mycolor1,'LineWidth',2,'HandleVisibility','off')
    xline(mean(diff_low(:),'omitmissing'),'--','Color',mycolor2,'LineWidth',2,'HandleVisibility','off')
    xlabel('error in width (m)')
    ylabel('#')
    legend('show')
    grid on;
    title(sprintf('%s original width of %0.1f %s %0.1f m',EARTH_RIVERS{i},mean(true_width(:),'omitmissing'),char(177),std(true_width(:),'omitmissing')))
   


end

figure;
hold on

for i = 1:3

    true_width = all_width{i,1}(:);
    high_width = all_width{i,2}(:);
    low_width  = all_width{i,3}(:);

    
    valid_high = ~isnan(true_width) & ~isnan(high_width);
    valid_low  = ~isnan(true_width) & ~isnan(low_width);
    %Background observations
    scatter(true_width(valid_high), high_width(valid_high), ...
        18, mycolor1, 'filled', ...
        'MarkerFaceAlpha',0.1, ...
        'MarkerEdgeAlpha',0.1, ...
        'HandleVisibility','off');

    scatter(true_width(valid_low), low_width(valid_low), ...
        18, mycolor2, 'filled', ...
        'MarkerFaceAlpha',0.1, ...
        'MarkerEdgeAlpha',0.1, ...
        'HandleVisibility','off');

    %Means
    x = mean(true_width,'omitnan');

    y_high = mean(high_width,'omitnan');
    y_low  = mean(low_width,'omitnan');

    sd_high = std(high_width,'omitnan');
    sd_low  = std(low_width,'omitnan');

    %175 m/pixel
    p = errorbar(x, y_high, sd_high, ...
        'o', ...
        'Color',mycolor1, ...
        'MarkerFaceColor',mycolor1, ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',10, ...
        'LineWidth',2.5, ...
        'CapSize',14, ...
        'DisplayName','175 m/pixel');

    %30 m/pixel
    p2 = errorbar(x, y_low, sd_low, ...
        'o', ...
        'Color',mycolor2, ...
        'MarkerFaceColor',mycolor2, ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',10, ...
        'LineWidth',2.5, ...
        'CapSize',14, ...
        'DisplayName','30 m/pixel');
end


plot([0 2100],[0 2100], ...
    '--', ...
    'Color','k', ...
    'LineWidth',1.5, ...
    'HandleVisibility','off');

%% Formatting
xlim([0 1500])
ylim([0 2100])
%axis equal
grid on
xlabel('True width (m)','FontSize',14)
ylabel('Measured width (m)','FontSize',14)



legend([p2(1) p(1)],'Location','northwest')

set(findall(gcf,'-property','FontSize'),'FontSize',12)