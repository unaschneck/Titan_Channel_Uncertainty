function [profile_brightness_all,profile_brightness_inside, profile_brightness_outside,contrast_all] = calc_brightness_profile_crossstream(x_center,y_center,x_left,y_left,x_right,y_right,sar_data,width,num_subdivisions)
% calculates the crossstream SAR brightness of the feature downstream along discrete cross-sections , including contrast with surronding landscapes
% INPUT:
%   x_center                   = x coordinates of the centerline of the feature in East-West dimensions (m)
%   y_center                   = y coordinates of the centerline of the feature in North-South dimensions (m)
%   x_left                     = x coordinates of the left bank of the feature in East-West dimensions (m)
%   y_left                     = y coordinates of the left bank of the feature in East-West dimensions (m)
%   x_right                    = x coordinates of the right bank of the feature in East-West dimensions (m)
%   y_right                    = y coordinates of the right bank of the feature in East-West dimensions (m)
%   sar_data                   = SAR brightness of the pass (2D array)
%   width                      = 1D array of width values along the centerline at each (x_center,y_center) (m)
%   num_subdivisions           = number of subdivisions in downstream direction to cut river into for a sense of the downstream variability in cross-stream brightness
% OUTPUT:
%   profile_brightness_all     = the average brightness of the feature, including its surrondings (3.5x local width)
%   profile_brightness_inside  = the average brightness of only pixels within the left and right banks
%   profile_brightness_outside = the average brightness outside of the pixels contained within the left and right bank (aka the surronding terrain)
%   contrast_all               = the average contrast between the brightness inside the feature and the average brightness in the nearby surronding terrain 


make_plot = 1;

n_downstream = 5; % number of downstream samples
n_crossstream = 20; % number of crosstream samples


[data_rows,data_cols] = size(sar_data);

[xgrid,ygrid] = meshgrid(1:data_cols,1:data_rows);


downstream_length = round(numel(x_center)/num_subdivisions); % 100 points / 5 subdivisions = 20 points downstream (10 points upstream center and 10 points downstream of center)
c_idx = linspace(round(downstream_length/2),numel(x_center)-round(downstream_length/2),num_subdivisions); % center of rectangles 

% compute normal relative to center of rectangle
dx = gradient(x_center);
dy = gradient(y_center);
% normalize
mag_tan = vecnorm([dx; dy]);
dx = dx./mag_tan;
dy = dy./mag_tan;
% normal is just perpendicular to this
nx = -dy;
ny = dx;

dx = dx(:);

dx = [dx(1); dx];
dy = [dy(1); dy];
nx = [nx(1); nx];
ny = [ny(1); ny];


b_inside = cell(num_subdivisions,1);
b_outside = cell(num_subdivisions,1);
b_all = cell(num_subdivisions,1);


% figure;
% imagesc(sar_data)
% hold on
% plot(x_left,y_left,'-b')
% plot(x_right,y_right,'-b')

if make_plot
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
end

for i = 1:num_subdivisions


        % polygon of channel shape 
        X_RIVER = [x_left(:); flipud(x_right(:))];
        Y_RIVER = [y_left(:); flipud(y_right(:))];

%
%
%               upstream
%       P1 ---------------- P2
%        |                  |
%        |                  |
%        |      x0,y0       |
%        |                  |
%        |                  |
%       P4 ---------------- P3
%            downstream
%
% each point is a distance 3/2*width away from the center


        % make rectangular region of interest
        c = c_idx(i);
        Tx = dx(c);
        Ty = dy(c);
        Nx = nx(c);
        Ny = ny(c);

        rec_diag = ((3/2)*width(c));
        W = rec_diag*cosd(45);
        L = rec_diag*sind(45);

        xc = x_center(c);
        yc = y_center(c);
        
        P1_X = xc - W*Nx + L*Tx;
        P1_Y = yc - W*Ny + L*Ty;
        
        P2_X = xc + W*Nx + L*Tx;
        P2_Y = yc + W*Ny + L*Ty;
        
        P3_X = xc + W*Nx - L*Tx;
        P3_Y = yc + W*Ny - L*Ty;
        
        P4_X = xc - W*Nx - L*Tx;
        P4_Y = yc - W*Ny - L*Ty;

        X_REGION = [P1_X P2_X P3_X P4_X P1_X];
        Y_REGION = [P1_Y P2_Y P3_Y P4_Y P1_Y];
       
        % filter points within rectangular region
        inside_region = inpolygon(xgrid,ygrid,X_REGION,Y_REGION);
        % filter points within river
        inside_channel = inpolygon(xgrid,ygrid,X_RIVER,Y_RIVER);

        % subset of data within just the rectangular region
        RECT = sar_data;
        RECT(~inside_region) = NaN;
        % subset of data within the rectangular region AND within the channel
        ROI = sar_data;
        ROI(~(inside_region & inside_channel)) = NaN;

        t_down = linspace(0,1,n_downstream);
        t_cross = linspace(0,1,n_crossstream);

        % points on left edge (P1 -> P4)
        xL = P1_X + t_down.*(P4_X - P1_X);
        yL = P1_Y + t_down.*(P4_Y - P1_Y);
        
        % points on right edge P2 -> P4)
        xR = P2_X + t_down.*(P3_X - P2_X);
        yR = P2_Y + t_down.*(P3_Y - P2_Y);

        % points on upstream edge (P1 -> P2)
        xT = P1_X + t_cross.*(P2_X - P1_X);
        yT = P1_Y + t_cross.*(P2_Y - P1_Y);
            
        % points on downstream edge (P4 -> P3)
        xB = P4_X + t_cross.*(P3_X - P4_X);
        yB = P4_Y + t_cross.*(P3_Y - P4_Y);

       
       b_inside{i} = nan(n_downstream,n_crossstream);
       b_outside{i} = nan(n_downstream,n_crossstream);

        % extract brightness along slices
        for s = 1:n_downstream % dont do exact edges because pixels can erroneously fall out of region?
            [px{i}(s,:),py{i}(s,:),b{i}(s,:)] = improfile(sar_data,[xL(s) xR(s)],[yL(s) yR(s)],n_crossstream);
       

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
        I_surronding =mean_crossstream_outside_per_sub{i};

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
            plot(x_center(c_idx(i)), y_center(c_idx(i)),'or','MarkerFaceColor','r')

            for ii = 1:n_downstream
                plot([xL(ii) xR(ii)],[yL(ii) yR(ii)],'-k','LineWidth',0.5)

                
              
                plot(px{i}(ii,inside{ii}), py{i}(ii,inside{ii}), 'sm', 'MarkerFaceColor','m')
                
                
                plot(px{i}(ii,~inside{ii}), py{i}(ii,~inside{ii}),'sb', 'MarkerFaceColor','b')
                

            end

            for ii = 1:n_crossstream
                plot([xT(ii) xB(ii)],[yT(ii) yB(ii)],'-k','LineWidth',0.5)
            end
           


            title('INSIDE SUBSECTION')
            axis xy
            colormap('bone')
            xlim([min(x_right,[],"all","omitmissing") - 100, max(x_right,[],"all","omitmissing") + 100]);
            ylim([min(y_right,[],"all","omitmissing") - 100, max(y_right,[],"all","omitmissing") + 100]);

            
            h2 = subplot(3,1,2);

            imagesc(sar_data,'AlphaData',0.7)
            hold on
            imagesc(ROI,'AlphaData',~isnan(ROI));
            
            plot(X_REGION, Y_REGION,'-k')
            plot(x_center(c_idx),y_center(c_idx),'or')
            plot(x_center(c_idx(i)), y_center(c_idx(i)),'or','MarkerFaceColor','r')


            for ii = 1:n_downstream
                plot([xL(ii) xR(ii)],[yL(ii) yR(ii)],'-k','LineWidth',1,'Color',line_color(ii,:))
            end
            for ii= 1:n_crossstream
                plot([xT(ii) xB(ii)],[yT(ii) yB(ii)],'-k','LineWidth',0.5)
            end


            title('INSIDE CHANNEL AND SUBSECTION')
            axis xy
            colormap('bone')

            xlim([min(x_right,[],"all","omitmissing") - 100, max(x_right,[],"all","omitmissing") + 100]);
            ylim([min(y_right,[],"all","omitmissing") - 100, max(y_right,[],"all","omitmissing") + 100]);
    
            h3 = subplot(3,1,3);

            grid on;
            hold on

            
            for ii = 1:n_downstream
            
                brightness_across = b{i}(ii,:);
                scatter(linspace(0,1,n_crossstream),brightness_across,40,~inside{ii},'filled')
                plot(linspace(0,1,n_crossstream),brightness_across,'-o','Color',line_color(ii,:),'LineWidth',2)
               

            end
            xlabel('cross stream fraction')
            ylabel('brightness')

            drawnow

        %   if i == 1
        %       gif('BigBend.gif','DelayTime',1)
        %   else
        %       gif;
        %   end
        % 
        % 
        % end


       
    

end

profile_brightness_inside = mean(cell2mat(mean_crossstream_inside_per_sub(:)),1,'omitnan');

profile_brightness_outside = mean(cell2mat(mean_crossstream_outside_per_sub(:)),1,'omitnan');

profile_brightness_all = mean(reshape(cell2mat(b_all),[],n_crossstream),1,'omitnan');

contrast_all = mean(contrast_per_sub,'omitnan');





end


