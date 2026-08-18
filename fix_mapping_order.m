function [x_center,y_center,x_left,y_left,x_right,y_right] = fix_mapping_order(x_center,y_center,x_left,y_left,x_right,y_right)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that each line is mapped in the same direction downstream. Centerline is used as the reference direction. If either 
% of the banks is reversed relative to the centerline mapping then it is flipped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   x_center  : x coordinate centerline in original mapping direction downstream
%   y_center  : y coordinate centerline in original mapping direction downstream
%   x_left    : x coordinate left bank in original mapping direction downstream
%   y_left    : y coordinate left bank in original mapping direction downstream
%   x_right   : x coordinate right bank in original mapping direction downstream
%   y_right   : y coordinate right bank in original mapping direction downstream
% OUTPUT 
%   x_center  : x coordinate centerline in original mapping direction downstream (unchanged)
%   y_center  : y coordinate centerline in original mapping direction downstream (unchanged)
%   x_left    : x coordinate left bank in centerling mapping direction downstream (flipped if needed)
%   y_left    : y coordinate left bank in centerling mapping direction downstream (flipped if needed)
%   x_right   : x coordinate right bank in centerling mapping direction downstream (flipped if needed)
%   y_right   : y coordinate right bank in centerling mapping direction downstream (flipped if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the distance from the first point on the left/right bank to the centerline's first point should always be shorter than the 
% distance to the centerline's last point. Otherwise, flip it like its hot

first_center = [x_center(1),y_center(1)];
first_leftbank = [x_left(1),y_left(1)];
last_leftbank = [x_left(end),y_left(end)];
first_rightbank = [x_right(1),y_right(1)];
last_rightbank = [x_right(end),y_right(end)];

first_left2cent = sqrt(sum((first_leftbank - first_center).^2));
last_left2cent = sqrt(sum((last_leftbank - first_center).^2));

first_right2center =  sqrt(sum((first_rightbank - first_center).^2));
last_right2center = sqrt(sum((last_rightbank - first_center).^2));

if first_left2cent > last_left2cent
    x_left = flip(x_left);
    y_left = flip(y_left);
    disp('flipped order of left bank to match centerline order')
end

if first_right2center > last_right2center
    x_right = flip(x_right);
    y_right = flip(y_right);
    disp('fipped order of left bank to match centerline order')
end

end