function [avg_sin] = calc_average_sinuosity(feature_name,max_subdivisions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the average sinosity for a river (with information stored in file 
% 'feature_name') for increasing number of subdivisions up to 'max_subdivisions'.
% Sinuosity is defined such that 1 describes the sinuosity of a  straight line
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   feature_name     : name of .mat file containing the feature information (e.g., 'EquatorialEx_feature_info') [string]
%   max_subdivisions : maximum number of subdivisions of river to consider [1 x 1]
% OUTPUT:
%   avg_sin          : average sinuosity for river with different number of subdivisions [1 x max_subdivisions] 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_plot = 0; % if set to 1 will show a plot of the river centerline and the average sinusoity for increasing segment lengths
test_straight_line = 0; % if set to 1 will override river information provided and run a unit testing for a straight line at a 45 degree angle (should give a sinuosity of 1 for all segment lengths)

river = load(feature_name); % e.g., load('EquatorialEx_feature_info')

river_mean_width = mean(river.geometry.width,'omitmissing');

x_centerline = river.geometry.center_dimensionless(:,1); % made dimensionless by dividing by the river width in the original extraction
y_centerline = river.geometry.center_dimensionless(:,2); % made dimensionless by dividing by the river width in the original extraction


d_downstream = river.geometry.distance;


if test_straight_line

    feature_name = 'straight_line_test_case';
    river_mean_width = 1;
    
    x_centerline = linspace(0,100,100);
    y_centerline = x_centerline;
    
    
    d_segment = sqrt(diff(x_centerline).^2 + diff(y_centerline).^2);
    d_downstream = [0 cumsum(d_segment)];

end


if show_plot
    figure('Name',strrep(feature_name,'_',' '))
    hold on;
    axis equal;
    xlabel('x [dimensionless]')
    ylabel('y [dimensionless]')
    title([strrep(feature_name,'_',' '), ' centerline'])
    grid on

    plot(x_centerline,y_centerline,'.-b')
end

for n = 1:max_subdivisions
    

    subdivision = linspace(0,d_downstream(end),n+1);

    for s = 1:numel(subdivision)-1
    
        cut = d_downstream >= subdivision(s) & d_downstream <= subdivision(s+1);
    
        x_seg = x_centerline(cut);
        y_seg = y_centerline(cut);
    
        d_river_path = sum(sqrt(diff(x_seg).^2 + diff(y_seg).^2))*river_mean_width;
    
        d_straight_path = sqrt((x_seg(end)-x_seg(1))^2 + (y_seg(end)-y_seg(1))^2)*river_mean_width;
    
        sin_at_n(s) = d_river_path/d_straight_path;
    
    end

    avg_sin(n) = mean(sin_at_n,'omitmissing');

end

    
if show_plot

    sin_distance_scale = d_downstream(end)./(1:max_subdivisions); % segment length = total channel length / number of subdivisions

    figure('Name','Sinuosity')


    plot(sin_distance_scale./1000,avg_sin,'-r','LineWidth',2);
    yline(avg_sin(1),'--','global sinuosity') 
    grid on
    xlabel('Segment Length (km)')
    ylabel('Sinuosity')
    title(strrep(feature_name,'_',' '))

end

end

