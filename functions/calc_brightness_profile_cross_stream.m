function [profile_brightness_all,profile_brightness_inside, profile_brightness_outside,contrast_all] = calc_brightness_profile_cross_stream(feature_name,lat_center,lon_center,lat_left,lon_left,lat_right,lon_right,ellipsoid,sar_data,sar_file,width,start_row,start_col,num_subdivisions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the crossstream SAR brightness of the feature downstream along discrete cross-sections , including contrast with surronding landscapes
% INPUT:
%   feature_name               = unique name of mapped feature (e.g., BigBend)
%   lat_center                 = latitude (deg) of centerline
%   lon_center                 = longitude (deg) of centerline
%   lat_left                   = latitude (deg) of left bank
%   lon_left                   = longitude (deg) of left bank
%   lat_right                  = latitude (deg) of right bank
%   lon_right                  = longitude (deg) of right bank
%   ellipsoid                  = ellipsoid to use (e.g., ellipsoid = [2575000,0] for Titan)
%   sar_data                   = SAR brightness of the pass (2D array)
%   sar_file                   = filename containing sar data
%   width                      = 1D array of width values along the centerline at each (x_center,y_center) (m)
%   start_row                  = index of first row of subset of data with river
%   start_col                  = index first column of subset of data with river
%   num_subdivisions           = number of subdivisions in downstream direction to cut river into for a sense of the downstream variability in cross-stream brightness
% OUTPUT:
%   profile_brightness_all     = the average brightness of the feature, including its surrondings (3.5x local width)
%   profile_brightness_inside  = the average brightness of only pixels within the left and right banks
%   profile_brightness_outside = the average brightness outside of the pixels contained within the left and right bank (aka the surronding terrain)
%   contrast_all               = the average contrast between the brightness inside the feature and the average brightness in the nearby surronding terrain 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
sanity_check = 1;

% coordinates of the centerline of the feature in East-West dimensions (m)
[x_center,y_center,~] = deg2meters(lat_center,lon_center,ellipsoid,lat_center(1),lon_center(1));
% coordinates of the left bank of the feature in East-West dimensions (m)
[x_left,y_left,up_left] = deg2meters(lat_left,lon_left,ellipsoid,lat_center(1),lon_center(1));
% coordinates of the right bank of the feature in East-West dimensions (m)
[x_right,y_right,up_right] = deg2meters(lat_right,lon_right,ellipsoid,lat_center(1),lon_center(1));

% row,col index of pixel coordinate for centerline
[x_center_pix, y_center_pix] = transform_coordinates(lat_center, lon_center, sar_file);
% row,col index of pixel coordinate for left bank
[x_left_pix, y_left_pix] = transform_coordinates(lat_left, lon_left, sar_file);
% row,col index of pixel coordinate for right bank
[x_right_pix, y_right_pix] = transform_coordinates(lat_right, lon_right, sar_file);

% translate indices to indices in subset
[y_center_pix,x_center_pix] = index_for_subset(y_center_pix,x_center_pix,start_row,start_col);
[y_left_pix,x_left_pix] = index_for_subset(y_left_pix,x_left_pix,start_row,start_col);
[y_right_pix,x_right_pix] = index_for_subset(y_right_pix,x_right_pix,start_row,start_col);

make_plot = 1;

n_downstream =  10; % number of downstream samples
n_crossstream = 11; % number of crosstream samples

dis_away_from_bank = 3; % factor of width away from bank (e.g., 3 = 3 widths away from the bank on either side)

[data_rows,data_cols] = size(sar_data);
[col_grid,row_grid] = meshgrid(1:data_cols,1:data_rows);

c_idx = zeros(num_subdivisions,1); 
dx_sub = zeros(num_subdivisions,1); % x component of downstream tangent for each subdivision
dy_sub = zeros(num_subdivisions,1); % y component of downstream tangent for each subdivision
nx_sub = zeros(num_subdivisions,1); % x component of downstream normal for each subdivision
ny_sub = zeros(num_subdivisions,1); % y component of downstream normal for each subdivision

edges = round(linspace(1, numel(x_center)+1, num_subdivisions+1));

% calculate tangent and normal for each subdivision
for i = 1:num_subdivisions
    
    % indexes of subdivision
    idx_start = edges(i);
    idx_end = edges(i+1)-1;
    idx = idx_start:idx_end;
    c_idx(i) = round(mean([idx_start idx_end]));

    % centerline coordinates for this subdivision
    x_local = x_center(idx);
    y_local = y_center(idx);

    % Remove the mean to get local orientation
    x_local = x_local - mean(x_local);
    y_local = y_local - mean(y_local);

    % first principal component gives representative tangent 
    % for the entire  subdivision rather than a tangent calculated 
    % at every point (slightly better than fitting line to endpoints
    % because less sensitive to jitter but aprox same most of the time
    [~,~,V] = svd([x_local(:), y_local(:)],'econ');
    tangent = V(:,1);
    % normalize
    tangent = tangent ./ norm(tangent);
    dx_sub(i) = tangent(1);
    dy_sub(i) = tangent(2);

    % make sure tangent points downstream
    downstream_direction = [x_center(idx_end) - x_center(idx_start), y_center(idx_end) - y_center(idx_start)];
    if dot([dx_sub(i),dy_sub(i)],downstream_direction) < 0
        % Reverse tangent if it points upstream
        dx_sub(i) = -dx_sub(i);
        dy_sub(i) = -dy_sub(i);
    end

    % compute the cross stream normal
    nx_sub(i) = -dy_sub(i);
    ny_sub(i) =  dx_sub(i);

end

if make_plot
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
end


for i = 1:num_subdivisions 


        % polygon of channel shape in pixel indices
        X_RIVER = [x_left_pix(:); flipud(x_right_pix(:))];
        Y_RIVER = [y_left_pix(:); flipud(y_right_pix(:))];

        % make rectangular region of interest (in meters and then convert back to pixel coordinates)
        c  = c_idx(i);
        Tx = dx_sub(i);
        Ty = dy_sub(i);
        Nx = nx_sub(i);
        Ny = ny_sub(i);
        
        xl = x_left(c);
        yl = y_left(c);
        xr = x_right(c);
        yr = y_right(c);
        xc = x_center(c);
        yc = y_center(c);

        W = dis_away_from_bank*(mean(width(:)));
        %W = dis_away_from_bank*max(width(:));%,"all","omitmissing");
 
        % determine which direction is OUTSIDE the river
       %  N = [Nx, Ny];
       %  left_bank_dir =   [xl - xc, yl - yc];
       %  right_bank_dir =  [xr - xc, yr - yc];
       %  left_side_sign = sign(dot(left_bank_dir, N));
       %  right_side_sign = sign(dot(right_bank_dir, N));
       % 
       % % bank point pushed W further out along the normal
       % leftC_X = xl + W * left_side_sign * Nx; 
       % leftC_Y = yl + W * left_side_sign * Ny; 
       % 
       % rightC_X = xr + W * right_side_sign * Nx; 
       % rightC_Y = yr + W * right_side_sign * Ny;

        % local direction from left bank to right bank
        bank_vec = [xr-xl, yr-yl];
        N_bank = bank_vec / norm(bank_vec);
        
        % make sure N_bank points in the same direction as normal
        if dot(N_bank,[Nx Ny]) > 0
            N_bank = -N_bank;
        end
        
        Nx = N_bank(1);
        Ny = N_bank(2);

        % make sure the rectangle is being extended from left to left edge and not from left to right edge
        if pdist([xl,yl;xr,yr],"euclidean") > pdist([xl-Nx,yl-Ny;xr,yr],"euclidean") 
            Nx = - Nx;
            Ny = - Ny;
        end

        % extend outward from the banks
        leftC_X  = xl - W*Nx;
        leftC_Y  = yl - W*Ny;
        
        rightC_X = xr + W*Nx;
        rightC_Y = yr + W*Ny;

        % each extended center pushed another W up and downstream along the tangent to form the vertices
        P1_X = rightC_X + W*Tx;   P1_Y = rightC_Y + W*Ty;   % right, downstream
        P3_X = rightC_X - W*Tx;   P3_Y = rightC_Y - W*Ty;   % right, upstream
        P2_X = leftC_X  + W*Tx;   P2_Y = leftC_Y  + W*Ty;   % left, downstream
        P4_X = leftC_X  - W*Tx;   P4_Y = leftC_Y  - W*Ty;   % left, upstream

        % figure;
        % hold on;
        % axis equal;
        % 
        % 
        % % Plot all centerline and banks
        % plot(x_center,y_center,'y-','LineWidth',1.5);
        % plot(x_left,y_left,'r-','LineWidth',1);
        % plot(x_right,y_right,'r-','LineWidth',1);
        % 
        % % Plot local centerline point
        % plot(xc,yc,'go','MarkerFaceColor','g','MarkerSize',8);
        % 
        % % Plot local banks
        % plot(xl,yl,'bo','MarkerFaceColor','b','MarkerSize',8);
        % plot(xr,yr,'bo','MarkerFaceColor','b','MarkerSize',8);
        % 
        % % Plot bank-to-bank line
        % plot([xl xr],[yl yr],'b-','LineWidth',3);

        total_cross_distance = pdist([rightC_X,rightC_Y;leftC_X,leftC_Y],"euclidean");
        if sanity_check
            leftbank2rightbank = pdist([xl,yl;xr,yr],"euclidean");
            leftbank2leftedge  = pdist([xl,yl;leftC_X,leftC_Y],"euclidean");
            rightbank2rightedge = pdist([xr,yr;rightC_X,rightC_Y],"euclidean");
            leftbank2rightedge  = pdist([xl,yl;rightC_X,rightC_Y],"euclidean");
            rightbank2leftedge = pdist([xr,yr;leftC_X,leftC_Y],"euclidean");

            if abs(leftbank2leftedge-W) > 1e-8 || abs(rightbank2rightedge - W) > 1e-8
                error(['Rectangle is not 3*W away from banks:\n' ...
                    'Average width: %0.0f\n' ...
                    'W = 3*avg width = %0.0f\n' ...
                    'distance from left bank to left edge: %0.0f\n' ...
                    'distance from right bank to right edge: %0.0f\n'], ...
                    mean(width(:)), ...
                    W, ...
                    leftbank2leftedge, ...
                    rightbank2rightedge)
            else
                % disp('all good for distance of rectangle relative to closest bank')
            end
            if abs((leftbank2rightedge - leftbank2rightbank)-W) > 1e-8  || abs((rightbank2leftedge - leftbank2rightbank)-W) > 1e-8
               error(['bank to opposite edge is not W + local width\n' ...
                   'Local width: %0.0f\n' ...
                   'Average width: %0.0f\n' ...
                   'W = 3*avg width = %0.0f\n' ...
                   'distance from left bank to right edge: %0.0f\n' ...
                   'distance from right bank to left edge: %0.0f\n' ...
                   'leftbank error: %0.0f\n' ...
                   'rightbank error: %0.0f\n'],leftbank2rightbank, ...
                   mean(width(:)), ...
                   W, ...
                   leftbank2rightedge, ...
                   rightbank2leftedge, ...
                   abs((leftbank2rightedge - leftbank2rightbank)-W), ...
                   abs((rightbank2leftedge - leftbank2rightbank)-W))
            else
                % disp('all good for distance of rectangle relative to opposite bank')
            end

        end

       
        % plot(leftC_X,leftC_Y,'ok');
        % plot(rightC_X,rightC_Y,'ok');
        % plot([xl leftC_X],[yl leftC_Y],'r-','LineWidth',3);
        % Plot local normal
        % quiver(xc,yc,Nx*W,Ny*W,0,'c','LineWidth',2);
        % 
        % % Plot tangent
        % quiver(xc,yc,Tx*W,Ty*W,0,'m','LineWidth',2);
        % 
        % % Plot rectangle
        % plot([P1_X P2_X P4_X P3_X P1_X], ...
        %      [P1_Y P2_Y P4_Y P3_Y P1_Y], ...
        %      'k-','LineWidth',3);
        % 
        % % Plot extended bank points
        % plot(leftC_X,leftC_Y,'ms','MarkerFaceColor','m','MarkerSize',8);
        % plot(rightC_X,rightC_Y,'ms','MarkerFaceColor','m','MarkerSize',8);


       % convert points in meter back to degrees and then to associated pixel index
       [P1_lat,P1_lon] = meters2deg(P1_X,P1_Y,up_right(c),ellipsoid,lat_center(1),lon_center(1));
       [P1_X_pix,P1_Y_pix] = transform_coordinates(P1_lat,P1_lon,sar_file);
        % convert points in meter back to degrees and then to associated pixel index
       [P2_lat,P2_lon] = meters2deg(P2_X,P2_Y,up_left(c),ellipsoid,lat_center(1),lon_center(1));
       [P2_X_pix,P2_Y_pix] = transform_coordinates(P2_lat,P2_lon,sar_file);
        % convert points in meter back to degrees and then to associated pixel index
       [P3_lat,P3_lon] = meters2deg(P3_X,P3_Y,up_right(c),ellipsoid,lat_center(1),lon_center(1));
       [P3_X_pix,P3_Y_pix] = transform_coordinates(P3_lat,P3_lon,sar_file);
        % convert points in meter back to degrees and then to associated pixel index
       [P4_lat,P4_lon] = meters2deg(P4_X,P4_Y,up_left(c),ellipsoid,lat_center(1),lon_center(1));
       [P4_X_pix,P4_Y_pix] = transform_coordinates(P4_lat,P4_lon,sar_file);

       [P1_Y_pix,P1_X_pix] = index_for_subset(P1_Y_pix,P1_X_pix,start_row,start_col);
       [P2_Y_pix,P2_X_pix] = index_for_subset(P2_Y_pix,P2_X_pix,start_row,start_col);
       [P3_Y_pix,P3_X_pix] = index_for_subset(P3_Y_pix,P3_X_pix,start_row,start_col);
       [P4_Y_pix,P4_X_pix] = index_for_subset(P4_Y_pix,P4_X_pix,start_row,start_col);

       X_REGION = [P1_X_pix P2_X_pix P4_X_pix P3_X_pix P1_X_pix];   % right-down -> left-down -> left-up -> right-up -> close rectangle
       Y_REGION = [P1_Y_pix P2_Y_pix P4_Y_pix P3_Y_pix P1_Y_pix];

       X_REGION = X_REGION(:);
       Y_REGION = Y_REGION(:);

       %       river (upstream)
       %       |  |
       %   P4--|--|---P3  rectangular region to sample over 
       %   |   |  |   |
       %   |   |  |   |
       %   P2--|--|---P1
       %       |  |
       %       river (downstream)


        % filter points within rectangular region
        inside_region = inpolygon(col_grid,row_grid,X_REGION,Y_REGION);
        % filter points within river
        inside_channel = inpolygon(col_grid,row_grid,X_RIVER,Y_RIVER);

        % subset of data within just the rectangular region
        RECT = sar_data;
        RECT(~inside_region) = NaN;
        % subset of data within the rectangular region AND within the channel
        ROI = sar_data;
        ROI(~(inside_region & inside_channel)) = NaN;

        t_down = linspace(0,1,n_downstream); 
        t_cross =linspace(0,1,n_crossstream);

        xperm_R = P1_X_pix + t_down.*(P3_X_pix - P1_X_pix);   % right side perimeter, down->up
        yperm_R = P1_Y_pix + t_down.*(P3_Y_pix - P1_Y_pix);

        xperm_L = P2_X_pix + t_down.*(P4_X_pix - P2_X_pix);   % left side, down->up
        yperm_L = P2_Y_pix + t_down.*(P4_Y_pix - P2_Y_pix);

        % points on upstream edge (P1 -> P2)
        xperm_T = P1_X_pix + t_cross.*(P2_X_pix - P1_X_pix);
        yperm_T = P1_Y_pix + t_cross.*(P2_Y_pix - P1_Y_pix);

        % points on downstream edge (P4 -> P3)
        xperm_B = P3_X_pix + t_cross.*(P4_X_pix - P3_X_pix);
        yperm_B = P3_Y_pix + t_cross.*(P4_Y_pix - P3_Y_pix);

        b_inside{i} = nan(n_downstream,n_crossstream);
        b_outside{i} = nan(n_downstream,n_crossstream);

        % extract brightness along slices
        for s = 1:n_downstream % dont do exact edges because pixels can erroneously fall out of region?
            [px{i}(s,:),py{i}(s,:),b{i}(s,:)] = improfile(sar_data,[xperm_L(s) xperm_R(s)],[yperm_L(s) yperm_R(s)],n_crossstream);


            inside_region = inpolygon(px{i}(s,:),py{i}(s,:),X_REGION,Y_REGION);
            inside_channel = inpolygon(px{i}(s,:),py{i}(s,:),X_RIVER,Y_RIVER);

            inside{s} = inside_region & inside_channel;


            b_inside{i}(s,inside{s}) = b{i}(s,inside{s});
            b_outside{i}(s,~inside{s}) = b{i}(s,~inside{s});

        end

        b_all{i} = b{i};

        % mean across downstream direction for this subsection
        mean_crossstream_inside_per_sub{i} = mean(b_inside{i},1,'omitnan');
        mean_crossstream_outside_per_sub{i} = mean(b_outside{i},1,'omitnan');
        mean_crosstrream_all{i} = mean(b{i});

        I_feature = mean_crossstream_inside_per_sub{i};
        I_surronding = mean_crossstream_outside_per_sub{i};

        contrast_per_sub(i) = calculate_contrast(mean_crossstream_inside_per_sub{i},mean_crossstream_outside_per_sub{i});

        if make_plot
            clf

            line_color = jet(n_downstream);

            h1 = subplot(3,1,1);
            imagesc(sar_data,'AlphaData',0.7);
            hold on;
            imagesc(RECT,'AlphaData',~isnan(RECT));

            plot(X_REGION, Y_REGION,'-k')

            plot(x_center(c_idx),y_center(c_idx),'or')
            plot(x_center(c), y_center(c),'or','MarkerFaceColor','r')

            for ii = 1:n_downstream
                % plot([xperm_L(ii) xperm_R(ii)],[yperm_L(ii) yperm_R(ii)],'-k','LineWidth',0.5)
                plot(px{i}(ii,inside{ii}), py{i}(ii,inside{ii}), 'sm', 'MarkerFaceColor','m')
                plot(px{i}(ii,~inside{ii}), py{i}(ii,~inside{ii}),'sb', 'MarkerFaceColor','b')
            end

            % for ii = 1:n_crossstream
            %     plot([xperm_T(ii) xperm_B(ii)],[yperm_T(ii) yperm_B(ii)],'-k','LineWidth',0.5)
            % end
            % plot(x_left_pix,y_left_pix,'.r')
            % plot(x_right_pix,y_right_pix,'.r')

            title('INSIDE SUBSECTION')
            axis xy
            colormap('bone')
            xlim([min(x_center_pix,[],"all","omitmissing") - 100, max(x_center_pix,[],"all","omitmissing") + 100]);
            ylim([min(y_center_pix,[],"all","omitmissing") - 100, max(y_center_pix,[],"all","omitmissing") + 100]);


            h2 = subplot(3,1,2);

            imagesc(sar_data,'AlphaData',0.7)
            hold on
            imagesc(ROI,'AlphaData',~isnan(ROI));

            plot(X_REGION, Y_REGION,'-k')
            plot(x_center(c_idx),y_center(c_idx),'or')
            plot(x_center(c_idx(i)), y_center(c_idx(i)),'or','MarkerFaceColor','r')


            for ii = 1:n_downstream
                plot([xperm_L(ii) xperm_R(ii)],[yperm_L(ii) yperm_R(ii)],'-k','LineWidth',1,'Color',line_color(ii,:))
            end
            for ii= 1:n_crossstream
                plot([xperm_T(ii) xperm_B(ii)],[yperm_T(ii) yperm_B(ii)],'-k','LineWidth',0.5)
            end

            % plot(x_left_pix,y_left_pix,'.r')
            % plot(x_right_pix,y_right_pix,'.r')

            title('INSIDE CHANNEL AND SUBSECTION')
            axis xy
            colormap('bone')

            xlim([min(x_center_pix,[],"all","omitmissing") - 100, max(x_center_pix,[],"all","omitmissing") + 100]);
            ylim([min(y_center_pix,[],"all","omitmissing") - 100, max(y_center_pix,[],"all","omitmissing") + 100]);

            h3 = subplot(3,1,3);

            grid on;
            hold on
            
           x_cross = linspace(-1,1,n_crossstream);

            for ii = 1:n_downstream
                
                % relative to midpoint of the points within the channel
                x_centered_m = x_cross.*total_cross_distance;

                brightness_across = b{i}(ii,:);

                scatter(x_centered_m./1000, brightness_across, 40, ~inside{ii}, 'filled')
            
                plot(x_centered_m./1000, brightness_across, '-o', 'Color', line_color(ii,:), 'LineWidth', 2)
            
            end

            % xline(0,'--k','Channel center','LineWidth',1.5)
            
            xlabel('cross stream distance from channel center (km)')
            ylabel('brightness')
            xlim([min(x_centered_m)/1000 max(x_centered_m)/1000])

            drawnow;


            if make_plot

                if i == 1
                    extract_flyby = extractBefore(sar_file, '_');
                    gif_name = fullfile(pwd,'results', [feature_name, '_',extract_flyby,'.gif']);
                    
                    gif(gif_name, 'DelayTime', 1);
                else
                    gif;
                end

            end

        end

end

profile_brightness_inside = mean(cell2mat(mean_crossstream_inside_per_sub(:)),1,'omitnan');
profile_brightness_outside = mean(cell2mat(mean_crossstream_outside_per_sub(:)),1,'omitnan');
profile_brightness_all = mean(reshape(cell2mat(b_all),[],n_crossstream),1,'omitnan');
contrast_all = mean(contrast_per_sub,'omitnan');


end

