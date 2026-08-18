function  [downstream_brightness_distance,brightness_profile]= calc_brightness_profile_downstream(x_feature,y_feature,sar_data,brightness_vals,x_left,y_left,x_right,y_right,downstream_distance,num_subdivisions)
%   calc along-channel brightness profile by measuring the average radar brightness within each subsection of the feature and
%   associating that measurement with the distance along the central path to the center of the subsection.
% INPUT:
%   x_feature                      : x coordinates of extracted SAR pixels (2D array of size brightness vals)
%   y_feature                      : y coordinates of extracted SAR pixels (2D array of size brightness vals)
%   brightness_vals                : SAR brightness of each extracted pixel (2D array)
%   x_left,y_left                  : left feature boundary coordinates
%   x_right,y_right                : right feature boundary coordinates
%   downstream_distance            : distance along feature centerline
%   num_subdivisions               : number of downstream sections
%
% OUTPUT:
%   downstream_brightness_distance : distance to center of each subsection downstream
%   brightness_profile             : mean brightness in each subsection

make_plot = 1;

downstream_brightness_distance = NaN(num_subdivisions,1);
brightness_profile = NaN(num_subdivisions,1);

downstream = round(linspace(1, numel(downstream_distance), num_subdivisions+1)); % divide downstream distance into equal sections

% Compute brightness for each downstream subsection
for i = 1:num_subdivisions

    if make_plot
        f1 = figure;
        hold on
        colormap bone
        axis equal
        axis xy
        
        
        plot(x_left,y_left,'.r','MarkerSize',10)
        plot(x_right,y_right,'.r','MarkerSize',10)
    end

    my_idx = downstream(i):downstream(i+1); % current subsection index along feature boundary

    % polygon bounding subsection
    x_box = [x_left(my_idx)' fliplr(x_right(my_idx)')];
    y_box = [y_left(my_idx)' fliplr(y_right(my_idx)')];
   
    inside_subsection = inpolygon(x_feature, y_feature, x_box,y_box);  % find all pixels inside subsection

    mid_point = round((downstream(i)+downstream(i+1))/2); % distance associated with subsection center

    downstream_brightness_distance(i) = downstream_distance(mid_point);

    if make_plot
        % Active polygon
        plot(x_box, y_box, '.b','MarkerSize',10)
    
        scatter(x_feature(inside_subsection), y_feature(inside_subsection), 20,  brightness_vals(inside_subsection), 'filled');
        colorbar
        %scatter(x_feature(~inside_subsection),y_feature(~inside_subsection),10,'rx') % for visual effect but kills the memory

        xlim([min(x_feature,[],"all","omitmissing") - 100, max(x_feature,[],"all","omitmissing") + 100]);
        ylim([min(y_feature,[],"all","omitmissing") - 100, max(y_feature,[],"all","omitmissing") + 100]);
        clim([0 max(brightness_vals,[],'all','omitmissing')])
    end

    % calc average brightness of pixels in subsection
    if any(inside_subsection,'all')

        brightness_profile(i) = mean(brightness_vals(inside_subsection),'omitnan');
        

    else

        brightness_profile(i) = NaN;

    end

   %fprintf('For %i subsections:\n',num_subdivisions)
   %fprintf('Cuts between indices %i - %i\n %i',min(my_idx),max(my_idx));

    if make_plot
        title(sprintf('River cut into %i piece(s) (mean brightness of %i of %i = %0.2f)',num_subdivisions,i,num_subdivisions,brightness_profile(i)))
    end
    
if make_plot    
    close(f1)
end

end



end