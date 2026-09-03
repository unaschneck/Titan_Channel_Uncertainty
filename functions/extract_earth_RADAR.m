function [lon,lat,radar_vals,brightness_vals] = extract_earth_RADAR(fname,x_left,y_left,x_right,y_right)


    [radar_vals, R] = readgeoraster(fname);
    radar_vals = double(radar_vals);
    
    lon_plot = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2));
    lat_plot = linspace(R.LatitudeLimits(2), R.LatitudeLimits(1), R.RasterSize(1));
    
    cols = 1:R.RasterSize(2);
    rows = 1:R.RasterSize(1);
    [lon, ~] = intrinsicToGeographic(R, cols, ones(size(cols)));
    [~, lat] = intrinsicToGeographic(R, ones(size(rows)), rows);
    
    
    % remove bad values values
    radar_vals(radar_vals < -2e38) = NaN;
    
    brightness_vals = radar_vals;

    figure;
    imagesc(lon_plot,lat_plot,log10(brightness_vals),'AlphaData', ~isnan(brightness_vals))
    colormap(bone)
    cb = colorbar();
    ylabel(cb,'Mean brightness (dB)','Rotation',270)
    grid on
    set(gca,'YDir','normal')
    hold on
    plot(x_left,y_left,'-r')
    plot(x_right,y_right,'-r')
    
    % Create 2D coordinate grids
    [Lon, Lat] = meshgrid(lon, lat);
    
    % River polygon
    river_poly_lon = [x_left(:); flipud(x_right(:))];
    river_poly_lat = [y_left(:); flipud(y_right(:))];
    
    % Pixels inside polygon
    river_pix = inpolygon(Lon, Lat, river_poly_lon, river_poly_lat);
    
    brightness_vals(~river_pix) = NaN;

 

end