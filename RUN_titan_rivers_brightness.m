clc
clear
close all
% restoredefaultpath
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAINS SCRIPT TO RUN TITAN BRIGHTNESS RIVER ANALYSIS FROM FRANK'S ORIGINAL CODEBASE 
% created 6/2026 by U. Schneck (ugschneck@gmail.com)
% This script will loop through all rivers with provided shapefiles and associated SAR
% data to determine the width, cross-stream and downstream brightness, and sinuoisty
% of the channel features. 
% Quick look figures produced:
%   (1) downstream channel shape with associated width 
%           "rivername_channel_path.png"
%   (2) mean cross channel brightness within, outside, and within + outside
%           "rivername_crosschannel_brightness.png"
%   (3) downstream brightness variation
%           "rivername_downstream_brightness.png"
%   (4) gif(s) of cross channel brightness calculation for each subset of reach
%           "rivername_tXX.gif"
% Summary data of mean values is saved to : Titan_rivers.xlsx 
% Full information for individual rivers is saved to : rivername_geometry.mat
%   - rivername from mapping (river_name)
%   - sinuosity of river (sinuosity)
%   - river width downstream (river_width) 
%   - downstream distance vector (downstream_distance)
%   - sar data array(s) used for mapping (sar_data)
%   - downstream brightness for each pass (down_brightness)
%   - cross stream brightness inside and outside river for each pass (all_brightness)
%   - cross stream brightness just inside river for each pass (inside_brightness)
%   - cross stream brighness just outside river for each pass (outside_brightness)
%   - contrast indices for each pass (contrast_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More in-depth optional features can be viewed for each function by settting make_plot
% within the function itself to true. Beware that this will create a lot of figures
% if run for all rivers to completion so probably better to run and pause on 
% created optional figures and save as needed before closing or MATLAB will 
% probably yell at you.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes advantage of external functions (see folder functions/external) including
%   (1) bidr_proj_v2       : projects coordinates into the BIDR data coordinates -- original author A. Hayes
%   (2) pds_label_parse_v2 : scrapes metadata from .cub files                    -- original author A. Hayes
%   (3) read_isis          : read in data from .cub files for passes             -- original author A. Hayes
%   (4) read_sbdr          : read in metadata from SBDR files                    -- original author A. Hayes
%   (5) bidr_proj_inv_v2   : inverts coordinates back to on the globe            -- original author A. Hayes
%   (5) batlow.mat         : good diverging colormap for figures                 -- https://zenodo.org/records/8409685
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script expects two sources of data
%   (1) BIDR_SBDR                : contains sar data and metadata for specific flybys (.cub, .sbdr)
%   (2) TITAN_MAPPED_CHANNELS    : mapped channel data, originally performed in ArcGIs -- original author F. Verdi (2023)
%   (3) counters                 : .mat file containing Lat, Lon, PassID information   -- original author A. Hayes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The mapped channels should be organized as:
%   Titan Features / FeatureName 
%                       |-------- LeftBank
%                       |-------- Centerline
%                       |-------- RightBank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(fullfile('functions')))                                               % all functions for extract and analysis channel data
addpath(genpath('data'))                                                              % raw flyby datas and marked feature data

% where is the data with the shape files
Titan_rivers_dir = fullfile('data','TITAN_MAPPED_CHANNELS');                          % location of directory with all Titan mapped river shape files
Titan_river = dir(Titan_rivers_dir);                                                  % all sub folders containing specific river data

addpath(genpath(Titan_rivers_dir));

fprintf('There are %i river features in %s\n',numel(Titan_river)-2,Titan_rivers_dir); % minus 2 to ignore folders '.' and '..'

results_dir = 'results';
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data from ArcGIS files -> .mat files 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ellipsoid = make_ellipsoid('Titan');                                                 % simple sphere with radius of Titan for geographic projection
num_points = 1000;                                                                   % number of points for the final mapped channels (will be evenly spaced)
max_num_subdivisions = 10;                                                           % maximum number of subdivisions to consider for cross and along channel brightness measurement

% Table of rivers with important values will be saved and reported at end of analysis
column_names = ["Feature" "Center Lat (deg)" "Center Lon (deg)" "Mean Width (m)" "Length Mapped (km)" "Mean Brightness" "Mean Contrast"  "SAR Track" "Sinuosity"];
column_type  = ["string" "double" "double" "string" "double" "double" "double"  "string" "double"];
table_sz     = [numel(Titan_river)-2 numel(column_names)];
RIVER_TABLE  = table(Size=table_sz,VariableTypes=column_type,VariableNames=column_names);
RIVER_TABLE.("Mean Brightness") = NaN(table_sz(1),1);                                % initialize brightness to NaN

k_table = 1;                                                                         % counter for river features (only reports those with shapefiles AND sar data)

% preallocate for memory
pp_all                          = cell(numel(Titan_river)-2,5);
pp_in                           = cell(numel(Titan_river)-2,5);
pp_out                          = cell(numel(Titan_river)-2,5);
contrast                        = cell(numel(Titan_river)-2,5);
downstream_brightness_distance  = cell(numel(Titan_river)-2,5);
brightness_profile_downstream   = cell(numel(Titan_river)-2,5);
avg_brightness                  = NaN(1,numel(Titan_river)-2);

disp('==================================================================')
disp('CALCULATING CHANNEL BRIGHTNESS METRICS')
disp('==================================================================')
for r = 1:numel(Titan_river)                                                           % loop through all feature data aka mapped channels (r ~ 26)


    if Titan_river(r).isdir && ~ismember(Titan_river(r).name , {'.','..'})             % ignore folders '.' and '..'

        feature_name = Titan_river(r).name;

        disp('==================================================================')


        RIVER_TABLE{k_table,1} = {feature_name};
        geo_river.river_name   = {feature_name};

        fprintf('Feature: %s\n',feature_name)

        % extract centerline and bank data from shapefile
        [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = extract_from_shp(feature_name,'Titan');
        
        % skip if there is no mapped centerline
        if isempty(x_center)
            continue
        end

        % make sure the lines are mapped in the same orientation downstream
        [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = fix_mapping_order(x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge);

        % interpolate the centerline to be evenly spaced
        [lat_center,lon_center] = even_geodesic_spacing(y_center,x_center,num_points,ellipsoid);
        % interpolate the left bank
        [lat_left,lon_left]     = even_geodesic_spacing(y_leftedge,x_leftedge,num_points,ellipsoid);
        % interpolate the right bank
        [lat_right,lon_right]   = even_geodesic_spacing(y_rightedge,x_rightedge,num_points,ellipsoid);

        % calculate total sinuosity of channel (of centerline)
        sin_river = calc_sinuosity(lon_center,lat_center,ellipsoid);

        % calculate river width along channel
        river_width = calc_river_width(lat_left,lon_left,lat_right,lon_right,ellipsoid);

        avg_width = mean(river_width,'omitmissing');                       % mean width of river 
        std_width = std(river_width,'omitmissing');                        % std of river width

        % normalize coordinates of centerline to average river width 
        [xd_center,yd_center] = dimensionless_coordinates(lat_center,lon_center,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));
        % normalize coordinates of left bank to average river width
        [xd_left,yd_left]     = dimensionless_coordinates(lat_left,lon_left,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));
        % normalize coordinates of right bank to average river width
        [xd_right,yd_right]   = dimensionless_coordinates(lat_right,lon_right,mean(river_width,'omitmissing'),ellipsoid,lat_center(1),lon_center(1));

        % find width of channel at each point the centerline
        width_centerline = width_along_centerline(lon_center,lat_center,lon_left,lat_left,river_width,ellipsoid);

        % find cumulative distance downstream
        downstream_distance = cumulative_distance_downstream(lat_center,lon_center,ellipsoid);

        % SAVE TO TABLE
        RIVER_TABLE{k_table,2} = mean(y_center,'omitmissing');                      % mean lat of river
        RIVER_TABLE{k_table,3} = mean(x_center,'omitmissing');                      % mean lon of river
        RIVER_TABLE{k_table,4} = string(sprintf('%.2f ± %.2f', avg_width, std_width));
        RIVER_TABLE{k_table,5} = downstream_distance(end)/1000;
        RIVER_TABLE{k_table,9} = sin_river;

        % PRINT TO CONSOLE
        fprintf('Center (lat,lon) : (%.2f, %.2f) deg\n', mean(y_center,'omitmissing'),mean(x_center,'omitmissing'));
        fprintf('Mean width: %0.1f %s %0.2f m\n',avg_width,char(177),std_width)
        fprintf('Sinuosity: %0.2f\n',sin_river)
        fprintf('Total mapped distance downstream: %0.1f km\n',downstream_distance(end)/1000)

        geo_river.sinuosity  = sin_river;

        % PLOTS
        channel_fig = figure('Name','Channel path and downstream width');
        subplot(2,1,1)
        scatter(lon_center,lat_center,10,1:length(lat_center),"filled")
        colormap("winter")
        grid on;
        xlabel('lon (deg)')
        ylabel('lat (deg)')
        
        subplot(2,1,2)
        std_top = river_width + std_width;
        std_bot = river_width - std_width;
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

        geo_river.river_width         = river_width;
        geo_river.downstream_distance = downstream_distance;

       
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RADAR BRIGHTNESS ALONG AND ACROSS CHANNEL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % find which swath is associated with the feature based on central coordinates of feature
        [feature_swaths_all,m_per_pixel,segment_number_all] = find_swaths(lat_center,lon_center,ellipsoid);
        [feature_swaths,feature_seg] = remove_low_quality_SAR(feature_swaths_all,segment_number_all);                    % Remove low quality SAR passes (radar modes 3 and 11 are High Res SAR)

        fprintf('swath resolution of ~ %i m per pixel\n',round(m_per_pixel));

        if ~isempty(feature_swaths)
            RIVER_TABLE{k_table,8} = {strjoin("t" + string(feature_swaths_all), ", ")}; 

            for cass_pass = 1:numel(feature_swaths)                                    % instead of user input in original code, just loop through all relevant passes for full collection

               swath_number = feature_swaths(cass_pass);
               seg_number   = feature_seg(cass_pass);
               fprintf('...checking data from t%i...\n',swath_number);

               sar_name = strcat('t',num2str(swath_number));

               [sar_file,sar_data] = load_SAR_data(sar_name,seg_number);


               if isempty(sar_file)                                                   % if there is no good data file in directory to use then skip this
                    downstream_brightness_distance{k_table,cass_pass} =  NaN(max_num_subdivisions,1);
                    brightness_profile_downstream{k_table,cass_pass}  =  NaN(max_num_subdivisions,1);

                    pp_all{k_table,cass_pass}   = NaN;
                    pp_in{k_table,cass_pass}    = NaN;
                    pp_out{k_table,cass_pass}   = NaN;
                    contrast{k_table,cass_pass} = NaN;
                   continue
               end

               % transform into meters relative to starting position in local east, north, up coordinates (ENU)
               [x_center, y_center]           = deg2meters(lat_center,lon_center,ellipsoid,lat_center(1),lon_center(1));
               [x_left, y_left]               = deg2meters(lat_left,lon_left, ellipsoid,lat_center(1),lon_center(1));
               [x_right, y_right]             = deg2meters(lat_right,lon_right,ellipsoid,lat_center(1),lon_center(1));

               % get pixel coordinates of paths in sar image
               [x_center_pix, y_center_pix]   = transform_coordinates(lat_center, lon_center, sar_file);
               [x_left_pix, y_left_pix]       = transform_coordinates(lat_left, lon_left, sar_file);
               [x_right_pix, y_right_pix]     = transform_coordinates(lat_right, lon_right, sar_file);

               % extract pixel brightness of channel feature
               [x_feature,y_feature,start_col,start_row,sar_data,brightness_vals] = extract_brightness(x_left_pix,y_left_pix,x_right_pix,y_right_pix,sar_data,sar_file);

               geo_river.sar_data{cass_pass} = sar_data;

               my_color = winter(max_num_subdivisions);

               fig_bright_down = figure('Name','Downstream brightness');
               ax_bright_down = axes(fig_bright_down);
               hold on;
               xlabel('downstream distance')
               ylabel('brightness')
               grid on;
               title(feature_name)

               fig_bright_across = figure('Name','Crosstream brightness');
               ax_bright_across = axes(fig_bright_across);
               hold on;
               xlabel('cross stream distance')
               ylabel('brightness')
               grid on;
               title(feature_name)

                if ~all(isnan(brightness_vals(:))) % check if this is a mismatch with location or if data has a bunch of NaNs

                    % compute along channel and across channel SAR brightness profile by dividing feature into downstream subsections and calculating the average radar brightness within each subsecitons and with surronding terrain
  
                    % DOWNSTREAM
                    [downstream_brightness_distance{k_table,cass_pass},brightness_profile_downstream{k_table,cass_pass}]= calc_brightness_profile_downstream(x_feature,y_feature,sar_data,brightness_vals,x_left_pix,y_left_pix,x_right_pix,y_right_pix,downstream_distance,max_num_subdivisions);
                    plot(ax_bright_down,downstream_brightness_distance{k_table,cass_pass},brightness_profile_downstream{k_table,cass_pass},'-o','MarkerFaceColor',my_color(1,:),'Color',my_color(1,:),'LineWidth',3,'DisplayName',sprintf('%i segments',max_num_subdivisions))
    
                    % CROSS STREAM
                    [pp_all{k_table,cass_pass},pp_in{k_table,cass_pass},pp_out{k_table,cass_pass},contrast{k_table,cass_pass}] = calc_brightness_profile_cross_stream(feature_name,lat_center,lon_center,lat_left,lon_left,lat_right,lon_right,ellipsoid,sar_data,sar_file,width_centerline,start_row,start_col,max_num_subdivisions);
                    plot(ax_bright_across,pp_all{k_table,cass_pass},'-k','LineWidth',3,'DisplayName','all')
                    plot(ax_bright_across,pp_in{k_table,cass_pass},'-b','LineWidth',3,'DisplayName','inside')
                    plot(ax_bright_across,pp_out{k_table,cass_pass},'-r','LineWidth',3,'DisplayName','outside')

                else
                    downstream_brightness_distance{k_table,cass_pass} =  NaN(max_num_subdivisions,1);
                    brightness_profile_downstream{k_table,cass_pass} =  NaN(max_num_subdivisions,1);
                    pp_all{k_table,cass_pass} = NaN;
                    pp_in{k_table,cass_pass} = NaN;
                    pp_out{k_table,cass_pass} = NaN;
                    contrast{k_table,cass_pass} = NaN;
                end
                hold off;

            end

            % report average brightness over feature for the table for all passses
            avg_b = brightness_profile_downstream(k_table, 1:size(brightness_profile_downstream,2));    % extract cells
            avg_b = avg_b(~cellfun('isempty', avg_b));                                                  % remove empty cells
            avg_b = [avg_b{:}];                                                                          % concatenate 

            avg_brightness(k_table) = mean(avg_b,'all' ,'omitnan');                                     % average over relevant swaths

            RIVER_TABLE{k_table,6} = avg_brightness(k_table);                                           % average over all relevant passes
            fprintf('Average pixel brightness: %0.1f dB\n',avg_brightness(k_table))

            contrast_avg = mean(cell2mat(contrast(k_table,:)),'all','omitmissing');

            RIVER_TABLE{k_table,7} = contrast_avg;
            fprintf('Average contrast with surrondings: %0.1f dB\n',contrast_avg)

            % saving to .mat files for each valid river
            geo_river.down_brightness    = avg_b;
            geo_river.down_brightness    = downstream_brightness_distance(k_table,:);
            geo_river.all_brightness     = pp_all(k_table,:);
            geo_river.inside_brightness  = pp_in(k_table,:);
            geo_river.outside_brightness = pp_out(k_table,:);
            geo_river.contrast_index     = contrast(k_table,:);

            k_table = k_table + 1;                                                                      % counter for table rows
        end

            % saving quick look figures
            saveas(channel_fig,fullfile(results_dir,[feature_name, '_channel_path.png']))
            saveas(fig_bright_down,fullfile(results_dir,[feature_name, '_downstream_brightness.png']))
            saveas(fig_bright_across,fullfile(results_dir,[feature_name, '_crosschannel_brightness.png']))
            save(fullfile(results_dir, [feature_name, '_geometry.mat']), 'geo_river');

    end

    close all % closes all figures for each river before moving on so a truly egregious amount of figures aren't created

end

writetable(RIVER_TABLE,fullfile('.',results_dir,'Titan_Rivers.xlsx'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES FOR GLOBAL MEAN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) MAKE A MAP OF BRIGHTNESS ON THE GLOBE
figure('Name','Channel brightness across globe');
scatter(RIVER_TABLE.("Center Lon (deg)")(:), RIVER_TABLE.("Center Lat (deg)")(:), 150, RIVER_TABLE.("Mean Brightness")(:), 'filled','MarkerEdgeColor','k')
ylabel(['Latitude (' char(176) ')'])
xlabel(['Longitude (' char(176) ')'])
colormap("copper")
cb = colorbar();
ylabel(cb,'Mean brightness (dB)','Rotation',270)
grid on
xlim([-180 180])
ylim([-90 90])
clim([0 1.5])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
fontname("arial")

% (2) MAKE A MAP OF CONTRAST ON THE GLOBE
figure('Name','Channel contrast across globe');
scatter(RIVER_TABLE.("Center Lon (deg)")(:), RIVER_TABLE.("Center Lat (deg)")(:), 150, RIVER_TABLE.("Mean Contrast")(:), 'filled','MarkerEdgeColor','k')
ylabel(['Latitude (' char(176) ')'])
xlabel(['Longitude (' char(176) ')'])
colormap(batlow)
cb = colorbar();
ylabel(cb,'Mean contrast','Rotation',270)
grid on
xlim([-180 180])
ylim([-90 90])
clim([-0.5 0.5])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
fontname("arial")