function dis_edge = calc_dis_edge(xc,yc,x_bank,y_bank,ellipsoid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating how far downstream the bank point is along the river centerline
% by finding the along-channel coordinate using a downstream nearest-neighbor
% search of segments (makes bank points 2D -> 1D). Instead of a bank point being something like
% (X,Y) it is now Z meters downstream on the centerline.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
%       xc        : x coordinate of centerline
%       yc        : y coordinate of centerline
%       x_bank    : x coordinate of bank (e.g., left bank)
%       y_bank    : y coordinate of bank (e.g., left bank)
%       ellipsoid : ellispoid of planet to calculate distances properly 
%   OUTPUT
%       dis_edge  : cumulative distance from the river bank point to the projected point on the centerline 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% segment lengths and cumulative downstream distance
seg_len = distance([yc(1:end-1) xc(1:end-1)], [yc(2:end)   xc(2:end)], ellipsoid);
seg_cumdis = [0; cumsum(seg_len)];

% change coordinates from lat lon to local centered at the mean
lat0 = mean(yc);
lon0 = mean(xc);
[xc_m,yc_m,~] = latlon2local( ...
    yc,xc, ...
    zeros(size(yc)), ...
    [lat0 lon0 0]);

[xb_m,yb_m,~] = latlon2local( ...
    y_bank,x_bank, ...
    zeros(size(y_bank)), ...
    [lat0 lon0 0]);

% calculate the segments and their size of the centerline
Ax  = xc_m(1:end-1);
Ay  = yc_m(1:end-1);
ABx = diff(xc_m);
ABy = diff(yc_m);
AB2 = ABx.^2 + ABy.^2;

% only check a few segments ahead and behind the bank to speed up computation
check_back    = 2;
check_forward = 8;

% prep for loop
dis_edge = zeros(numel(x_bank),1);
last_seg = 1;

% loop over all the bank points
for i = 1:numel(x_bank)

    Px = xb_m(i);
    Py = yb_m(i);

    % search only nearest segments
    segment_1 = max(1,last_seg-check_back);
    segment_2 = min(numel(xc)-1,last_seg+check_forward);

    best_dist_fin = inf;
    best_seg = last_seg;
    best_frac_pos = 0;

    for k = segment_1:segment_2

        APx = Px - Ax(k);
        APy = Py - Ay(k);

        % if frac_seg = 0, then the point is projected exactly onto pt A
        % if frac_seg = 1, then the point is projected exactly onto pt B
        % otherwise, frac_seg is between 0 and 1 aka between A and B
        frac_seg = (APx*ABx(k) + APy*ABy(k))/AB2(k); % fraction between A and B
        frac_seg = max(0,min(1,frac_seg));

        Qx = Ax(k) + frac_seg*ABx(k);
        Qy = Ay(k) + frac_seg*ABy(k);

        d2 = (Px-Qx)^2 + (Py-Qy)^2;

        if d2 < best_dist_fin
            best_dist_fin = d2;
            best_seg = k;
            best_frac_pos = frac_seg;
        end
    end

    dis_edge(i) = seg_cumdis(best_seg) + best_frac_pos*seg_len(best_seg);

    % starts next points search on current best segment
    last_seg = best_seg;

end



end