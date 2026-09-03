function [lat,lon] = bidr_proj_inv_v2(samples,lines,sar_file); %lat,lon is degrees

%get the header info
if ~isstruct(sar_file);
    pos = strfind(sar_file,'.raw');
    if length(pos) > 0;
        info = parse_raw(sar_file);
    else;
        info = pds_label_parse_v2(sar_file);
    end;
else;
    info = sar_file;
end;

% %define global constants
% R = 2575000; %Titan's radius in meters

%get the map projection section of header
if isfield(info,'qube');
    infom = info.qube;
elseif isfield(info,'isiscube');
    infom = info.isiscube;
elseif isfield(info,'image_map_projection');
    infom = info.image_map_projection;
elseif isfield(info,'israw');
    infom = info;
else;
    disp('Unsupported File (Not PDS, ISIS 2, ISIS 3, or RAW!');
    line = nan;
    sample = nan;
    return;
end;

%if planet input use default radius
if ~isfield(infom,'aaxisradius');    
    if exist('planet','var');
        switch planet
            case{'Titan'};
                infom.aaxisradius = 2575;
            case{'Mars'};
                infom.aaxisradius = 3396.190;
        end;        
    elseif isfield(infom,'equatorial_radius');
        infom.aaxisradius = 2575;        
    end;    
end;

%define global constants
if isfield(infom,'aaxisradius');
    R = infom.aaxisradius.*1000; %radius is meters
else;
    %default to Titan    
    R = 2575000; %Titan's radius in meters
end;

if isfield(infom,'res1');
    map_res1 = 1./(360./2./pi./R.*infom.res1);
    map_res2 = 1./(360./2./pi./R.*(-infom.res2));
    infom.mapresolution = (map_res1 + map_res2)./2;
elseif isfield(infom,'scale');
    infom.mapresolution = infom.scale;
end;

%get map projection type
if ~isfield(infom,'mapprojectiontype');
    infom.mapprojectiontype = infom.projectionname;
end;
if strcmp(infom.mapprojectiontype,'obliquecylindrical');
    infom.mapprojectiontype = 'oblique cylindrical';
end;
proj_type = infom.mapprojectiontype;
pos = strfind(proj_type,'_');
if length(pos) > 0; proj_type(pos) = ' '; end;

%get required info from header and define required constants
switch proj_type
    case{'oblique cylindrical'};
        if isfield(info,'isiscube');
            line_proj_offset = infom.upperleftcornery./infom.pixelresolution - 0.5;
            sample_proj_offset = -(infom.upperleftcornerx./infom.pixelresolution + 0.5);
        else
            line_proj_offset = infom.lineprojectionoffset;
            sample_proj_offset = infom.sampleprojectionoffset;
        end;
        if isfield(infom,'polelatitude');
            infom.obliqueprojpolelatitude = infom.polelatitude;
        end;
        if isfield(infom,'polelongitude');
            infom.obliqueprojpolelongitude = infom.polelongitude;
        end;
        if isfield(infom,'polerotation');
            infom.obliqueprojpolerotation = infom.polerotation;
        end;
        lat_p = infom.obliqueprojpolelatitude;
        lon_p = infom.obliqueprojpolelongitude;
        rot_p = infom.obliqueprojpolerotation;
        map_res = infom.mapresolution;       
        %create numbers from label value string
        vars = {'line_proj_offset','sample_proj_offset',...
            'lat_p','lon_p','rot_p','map_res'};
        for j=1:length(vars);
            if isstr(eval(vars{j}));
                eval(['pos = strfind( ' '' vars{j} '' ',''<'');']);
                if length(pos) > 0;
                    eval([vars{j} ' = str2num( ' ...
                        vars{j} '(1:(pos(1)-1)) );']);
                else;
                    eval([vars{j} ' = str2num( ' vars{j} ');']);
                end;
            end;
        end;
        %%fix the sign of the ISIS 2 cube offsets!!!
        if isfield(info,'qube');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
        %%PROJECTION TRANSFER EQUATIONS OBTAINED FROM PDS SIS
        %switch to pos east lon
        lon_p = -lon_p;
        if lon_p < -180; lon_p = lon_p +360; end;

        %convert to radians
        lon_p = lon_p*pi/180;
        lat_p = lat_p*pi/180;
        rot_p = rot_p*pi/180;
        
        %define M
        m1 = [ cos(rot_p) sin(rot_p) 0 ; ...
            -sin(rot_p) cos(rot_p) 0 ; ...
            0          0          1 ];
        m2 = [ cos(pi/2 - lat_p) 0 -sin(pi/2-lat_p) ; ...
            0                 1  0               ; ...
            sin(pi/2 - lat_p) 0  cos(pi/2-lat_p) ];
        m3 = [ cos(lon_p) sin(lon_p) 0 ; ...
            -sin(lon_p) cos(lon_p) 0 ; ...
            0          0          1 ];
        M = m1*m2*m3;
    case{'polarstereographic'};
        if isfield(info,'isiscube');
            line_proj_offset = infom.upperleftcornery./infom.pixelresolution - 0.5;
            sample_proj_offset = -(infom.upperleftcornerx./infom.pixelresolution + 0.5);
        else
            line_proj_offset = infom.lineprojectionoffset;
            sample_proj_offset = infom.sampleprojectionoffset;
        end;        
        center_lon = infom.centerlongitude;
        center_lat = infom.centerlatitude;
        if ~isfield(infom,'mapresolution');
            map_res = infom.scale;
        else;
            map_res = infom.mapresolution;
        end;
        if ~isfield(infom,'mapscale');
            map_scale = infom.pixelresolution./1000;
        else;
            map_scale = infom.mapscale;
        end;
        %create numbers from label value string
        vars = {'line_proj_offset','sample_proj_offset',...
            'center_lat','center_lon','map_res'};
        for j=1:length(vars);
            if isstr(eval(vars{j}));
                eval(['pos = strfind( ' '' vars{j} '' ',''<'');']);
                if length(pos) > 0;
                    eval([vars{j} ' = str2num( ' ...
                        vars{j} '(1:(pos(1)-1)) );']);
                else;
                    eval([vars{j} ' = str2num( ' vars{j} ');']);
                end;
            end;
        end;
        %convert to radians
        center_lat = center_lat.*pi/180;
        center_lon = center_lon.*pi/180;
        %%fix the sign of the PDS cube offsets!!!
        if isfield(info,'qube');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
        if abs(fix( sample_proj_offset ) - sample_proj_offset) == 0.5;
            soff = 0.5;
        else;
            soff = 1;
        end;
        loff = 1;
    case{'equirectangular cylindrical'};
        line_proj_offset = infom.lineprojectionoffset;
        sample_proj_offset = infom.sampleprojectionoffset;       
        center_lon = infom.centerlongitude;
        center_lat = infom.centerlatitude;
        map_res = infom.mapresolution;
        %create numbers from label value string
        vars = {'line_proj_offset','sample_proj_offset',...
            'center_lat','center_lon','map_res'};
        for j=1:length(vars);
            if isstr(eval(vars{j}));
                eval(['pos = strfind( ' '' vars{j} '' ',''<'');']);
                if length(pos) > 0;
                    eval([vars{j} ' = str2num( ' ...
                        vars{j} '(1:(pos(1)-1)) );']);
                else;
                    eval([vars{j} ' = str2num( ' vars{j} ');']);
                end;
            end;
        end;
        %convert to radians
        center_lat = center_lat.*pi/180;
        center_lon = center_lon.*pi/180;
        %%fix the sign of the PDS cube offsets!!!
        if isfield(info,'qube');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
         if isfield(info,'israw');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
        if abs(fix( sample_proj_offset ) - sample_proj_offset) == 0.5;
            soff = 0.5;
        else;
            soff = 1;
        end;
        loff = 1;
    case{'equirectangular'}
         disp('break');
        line_proj_offset = infom.upperleftcornery./infom.pixelresolution;
        sample_proj_offset = -infom.upperleftcornerx./infom.pixelresolution;       
        center_lon = infom.centerlongitude;
        center_lat = infom.centerlatitude;
        map_res = infom.mapresolution;
        %create numbers from label value string
        vars = {'line_proj_offset','sample_proj_offset',...
            'center_lat','center_lon','map_res'};
        for j=1:length(vars);
            if isstr(eval(vars{j}));
                eval(['pos = strfind( ' '' vars{j} '' ',''<'');']);
                if length(pos) > 0;
                    eval([vars{j} ' = str2num( ' ...
                        vars{j} '(1:(pos(1)-1)) );']);
                else;
                    eval([vars{j} ' = str2num( ' vars{j} ');']);
                end;
            end;
        end;
        %convert to radians
        center_lat = center_lat.*pi/180;
        center_lon = center_lon.*pi/180;
        %%fix the sign of the PDS cube offsets!!!
        if isfield(info,'qube');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
         if isfield(info,'israw');
            line_proj_offset = -line_proj_offset;
            sample_proj_offset = -sample_proj_offset;
        end;
        if abs(fix( sample_proj_offset ) - sample_proj_offset) == 0.5;
            soff = 0.5;
        else;
            soff = 1;
        end;
        loff = 1;
        
    otherwise
        disp(['Unsupported Map Projection: ' proj_type ' !!']);
        lat = nan;
        lon = nan;
        return
end;

%for i=1:length(lines);
    switch proj_type
        case{'oblique cylindrical'};
            %get the lat/lon
            line = lines;
            sample = samples;
            
            %derive lat_a and lon_a
            lon_a = (line - 1 - line_proj_offset)./map_res;
            lat_a = (sample - 1 - sample_proj_offset)./map_res;
            lon_a = lon_a.*pi/180;
            lat_a = lat_a.*pi/180;

            %derive lat and lon from lon_a and lat_a
            lat = asin( sin(lat_p).*sin(lat_a) - ...
                        cos(lat_p).*cos(lat_a).*cos(lon_a + rot_p) ).*180/pi;
            val = atan2( cos(lat_a).*sin(lon_a+rot_p),...
                sin(lat_p).*cos(lat_a).*cos(lon_a+rot_p)+cos(lat_p).*sin(lat_a));            
            lon = val.*180/pi + lon_p.*180/pi;
        case{'polarstereographic'}
            if ~isfield(infom,'positivelongitudedirection');
                if strcmp(infom.longitudedirection,'positiveeast');
                    infom.positivelongitudedirection = 'east';
                else;
                    infom.positive.ongitudedirection = 'west';
                end;
            end;
            if strcmp(infom.positivelongitudedirection,'east');                           
%                 lat_cen = pi./2;
%                 lon_cen = center_lon;
%                 x = (-1).*(sample_proj_offset - samples + soff)./map_res;
%                 y = (1).*(line_proj_offset - lines + loff)./map_res;
% %                 p = sqrt(x.^2 + y.^2);
% %                 c = 2.*atan2(p,2);
% %                 lat = asin( cos(c).*sin(lat_cen) + ...
% %                     (y.*sin(c).*cos(lat_cen))./p ).*180./pi;
% %                 lon = center_lon + atan2( (x.*sin(c)), ...
% %                     (p.*cos(lat_cen).*cos(c)-y.*sin(lat_cen).*sin(c))).*180./pi;
%                lon = atan2(x,-y)+lon_cen;
%               % lat = pi/2 - 2.*atan2(x,sin(lon-lon_cen));
%               % lat = pi/2 - 2.*atan2(-y,cos(lon-lon_cen));
%                lat = pi./2 - 2.*atan(sqrt(x.^2 + y.^2) );
%                lon = lon.*180./pi;
%                lat = lat.*180./pi;
                   mstruct = defaultm('stereo');
                   mstruct.origin = [center_lat.*180./pi ...
                       center_lon.*180./pi];
                   mstruct.geoid = [2575000 0];
                   mstruct = defaultm(mstruct);
%                    x = (samples-(2.0.*soff) - sample_proj_offset).*map_scale.*1000;
%                    y = (1).*(line_proj_offset - lines + (1.0.*loff)).*map_scale.*1000;
                    x = (samples-(1.0) - sample_proj_offset).*map_scale.*1000;
                    y = (1).*(line_proj_offset - lines + (1.0)).*map_scale.*1000;
                   [lat,lon] = minvtran(mstruct,x,y);                
            elseif strcmp(infom.positivelongitudedirection,'west');
%                 lon = atan2( (sample_proj_offset - samples + 1), ...
%                     (line_proj_offset - lines + 1) ).*180/pi + ...
%                     center_lon.*180/pi+center_lon.*180/pi./2;
%                 lat = atan2( (line_proj_offset - lines + 1),...
%                     (cos(lon.*pi/180).*map_res ) ).*180./pi - ...
%                     center_lat.*180/pi./2;
                disp(['West Longitude Not Supported for Stereographic Projection! '...
                    infom.positivelongitudedirection]);
                lat = nan;
                lon = nan;
                return
            else;
                disp(['Error in Positive Longitude Direction! '...
                    infom.positivelongitudedirection]);
                lat = nan;
                lon = nan;
                return
            end;
        case{'equirectangular'}
            sample_proj_offset = sample_proj_offset;
            if exist('map_res1','var');
                lat = (line_proj_offset - lines + loff)./map_res2;
            else;
                lat = (line_proj_offset - lines + loff)./map_res;
            end;
            if ~isfield(infom,'positivelongitudedirection');
                infom.positivelongitudedirection='east';
            end;
            if strcmp(infom.positivelongitudedirection,'east');
                if exist('map_res2','var');
                    lon = (-1).*(sample_proj_offset - samples +soff)...
                        ./map_res1./cos(center_lat) + center_lon.*180/pi;
                else;
                    lon = (-1).*(sample_proj_offset - samples +soff)...
                        ./map_res./cos(center_lat) + center_lon.*180/pi;
                end;
            elseif strcmp(infom.positivelongitudedirection,'west');
                if exist('map_res2','var');
                    lon = (sample_proj_offset - samples +soff)...
                        ./map_res2./cos(center_lat) + center_lon.*180/pi;
                else;
                    lon = (sample_proj_offset - samples +soff)...
                        ./map_res./cos(center_lat) + center_lon.*180/pi;
                end;
            else;
                disp(['Error in Positive Longitude Direction! '...
                    infom.positivelongitudedirection]);
                lat = nan;
                lon = nan;
                return
            end; 
        
        case{'equirectangular cylindrical'}
            sample_proj_offset = sample_proj_offset;
            if exist('map_res1','var');
                lat = (line_proj_offset - lines + loff)./map_res2;
            else;
                lat = (line_proj_offset - lines + loff)./map_res;
            end;
            if ~isfield(infom,'positivelongitudedirection');
                infom.positivelongitudedirection='east';
            end;
            if strcmp(infom.positivelongitudedirection,'east');
                if exist('map_res2','var');
                    lon = (-1).*(sample_proj_offset - samples +soff)...
                        ./map_res1./cos(center_lat) + center_lon.*180/pi;
                else;
                    lon = (-1).*(sample_proj_offset - samples +soff)...
                        ./map_res./cos(center_lat) + center_lon.*180/pi;
                end;
            elseif strcmp(infom.positivelongitudedirection,'west');
                if exist('map_res2','var');
                    lon = (sample_proj_offset - samples +soff)...
                        ./map_res2./cos(center_lat) + center_lon.*180/pi;
                else;
                    lon = (sample_proj_offset - samples +soff)...
                        ./map_res./cos(center_lat) + center_lon.*180/pi;
                end;
            else;
                disp(['Error in Positive Longitude Direction! '...
                    infom.positivelongitudedirection]);
                lat = nan;
                lon = nan;
                return
            end;            
        otherwise
            disp(['Unsupported Map Projection: ' proj_type ' !!']);
            lat = nan;
            lon = nan;
            return
    end;
%end;
