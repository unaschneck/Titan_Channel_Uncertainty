function [x_sar,y_sar,x_min,y_min,sub_sar,brightness_vals] = extract_brightness(x_left,y_left,x_right,y_right,sar_data,sar_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracts the brightness of the feature pixels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   x_left           : row pixel coordinates of left bank
%   y_left           : col pixel coordinate of left bank
%   x_right          : row pixel coordinates of right bank
%   y_right          : col pixel coordinates of right bank
%   sar_data         : normalized radar backscatter values 
%   sar_file         : file name with normalized radar backscatter values
% OUTPUTS
%   x_sar            : grid of x coordinates for brightness pixels in subset of sar data with a river
%   y_sar            : grid of y coordinates for brightness pixels in subset of sar data with a river
%   x_min            : starting column for subset of full data
%   y_min            : starting row for subset of full data
%   sub_sar          : subset grid of sar pixels containing river
%   brightness_vals  : brightness values of subset of grid containing river
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_plot = 1;

pad = 100; 

if x_left
% take a subset of data so don't hit memory issues
x_min = max(1, floor(min(x_left(:))) - pad);
x_max = min(size(sar_data,2), ceil(max(x_left(:))) + pad);
y_min = max(1, floor(min(y_left(:))) - pad);
y_max = min(size(sar_data,1), ceil(max(y_left(:))) + pad);

if y_min >= y_max || x_min >= x_max % river is not in scope of data? Shouldn't this be caught by counters?
    warning('river coordinates outside scope of %s\n',strrep(sar_file,'_',' '))
    x_sar = NaN;
    y_sar = NaN;
    x_min = NaN;
    y_min = NaN;
    sub_sar = NaN;
    brightness_vals = NaN;
    return
end

sub_sar = sar_data(y_min:y_max, x_min:x_max);


% Coordinates of each pixel in the ORIGINAL SAR coordinate system
[x_sar, y_sar] = meshgrid(x_min:x_max, y_min:y_max);

river_poly_x = [x_left(:); flipud(x_right(:))];
river_poly_y = [y_left(:); flipud(y_right(:))];

river_pix = inpolygon(x_sar,y_sar,river_poly_x,river_poly_y);


brightness_vals = sub_sar;
brightness_vals(~river_pix) = NaN;
    
    if make_plot

       

        figure('Name','Extracted Brightness')
        subplot(2,1,1)
        hold on
        imagesc(x_min:x_max, y_min:y_max, sub_sar)
        colorbar
        colormap bone
        axis equal
        axis xy
        plot(x_left,y_left,'.r','MarkerSize',2)
        plot(x_right,y_right,'.r','MarkerSize',2)
   
        
   
        xlim([x_min x_max]);
        ylim([y_min y_max]);

        
        colorbar   
    
        title(sprintf('Extracted Brightness Pixels (%s)',strrep(sar_file,'_',' ')))
    
        
        subplot(2,1,2)
        hold on
        imagesc(x_min:x_max, y_min:y_max, brightness_vals, 'AlphaData', ~isnan(brightness_vals))
        colorbar
        colormap bone
        axis equal
        axis xy
        xlim([x_min x_max]);
        ylim([y_min y_max]);


    end



end