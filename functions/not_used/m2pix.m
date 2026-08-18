function distance_pix = m2pix(distance_m,m_per_pixel)
% converts values in an array original in meters to values in number of pixels.
% e.g., width = [10 20 30] m for a resolution of 2 m/pixel would return
% a new array of widths = [5 4 6] pixels
% 
% INPUT:
%   distance_m   = array of distances in meters
%   m_per_pixel  = resolution of image (number of meters per one pixel)
% OUTPUT
%   distance_pix = array of distances in # of pixels

    distance_pix = round(distance_m./m_per_pixel);
    distance_pix(distance_pix==0) = 1;

end