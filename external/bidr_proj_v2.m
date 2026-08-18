function [line,sample,info,map_scale] = bidr_proj_v2(lats,lons,sar_file,proj_type,planet); %lat,lon is degrees

%convert lat/lon to radians
lats = lats.*pi/180;
lons = lons.*pi/180;

%get the header info
if ~isstruct(sar_file)
    pos = strfind(sar_file,'.raw');
    if length(pos) > 0
        info = parse_raw(sar_file);
    else
        info = pds_label_parse_v2(sar_file);
    end
else
    info = sar_file;
end

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
elseif isfield(infom,'equatorialradius');
    if isfield(infom,'equatorialradius_units');
        if strcmp(infom.equatorialradius_units,'meters');
            R = infom.equatorialradius;
        else;
            R = infom.equatorialradius.*1000;
        end;
    else;
        R = infom.equatorialradius.*1000;
    end;
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

%HACK THAT MAY OR MAY NOT WORK!!!!!
if strcmp(proj_type,'polarstereographic');
    proj_type = 'polar stereographic';
    infom.mapres = infom.scale;
    infom.mapscale = infom.pixelresolution./1000;
end;

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
    case{'polar stereographic'};
        if isfield(info,'isiscube');
            line_proj_offset = infom.upperleftcornery./infom.pixelresolution - 0.5;
            sample_proj_offset = -(infom.upperleftcornerx./infom.pixelresolution + 0.5);
        else
            line_proj_offset = infom.lineprojectionoffset;
            sample_proj_offset = infom.sampleprojectionoffset;
        end;        
        center_lon = infom.centerlongitude;
        center_lat = infom.centerlatitude;
        map_res = infom.mapresolution;
        map_scale = infom.mapscale;
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
    case{'equirectangular'};
        %disp('break');
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
    otherwise
        disp(['Unsupported Map Projection: ' proj_type ' !!']);
        lat = nan;
        lon = nan;
        return
end;
%disp(proj_type)
switch proj_type
    case{'oblique cylindrical'};
        %define XB        
        XB = zeros(3,numel(lats));
        XB(1,:) = R.*cos(lats).*cos(lons);
        XB(2,:) = R.*cos(lats).*sin(lons);
        XB(3,:) = R.*sin(lats);
        
        %Calculate XA
        XA = zeros(3,numel(lats));
        XA(1,:) = M(1,1).*XB(1,:)+M(1,2).*XB(2,:)+XB(3,:).*M(1,3);
        XA(2,:) = M(2,1).*XB(1,:)+M(2,2).*XB(2,:)+XB(3,:).*M(2,3);
        XA(3,:) = M(3,1).*XB(1,:)+M(3,2).*XB(2,:)+XB(3,:).*M(3,3);
        
        %calcualte lat_a and lon_a
        lon_a0 = atan2( XA(2,:),XA(1,:) ).*180./pi;
        lat_a0 = atan2( XA(3,:),sqrt( XA(1,:).^2 + XA(2,:).^2 ) ).*180./pi;
        
        %calculate alternative values for lon_a and lat_a
        lat_a = asin( sin(lat_p).*sin(lats) + ...
            cos(lat_p).*cos(lats).*cos(lons-lon_p) ).*180./pi;
        val = atan2(cos(lats).*sin(lons-lon_p), ...
            sin(lat_p).*cos(lats).*cos(lons-lon_p) - cos(lat_p).*sin(lats));
        lon_a = val.*180./pi - rot_p.*180./pi;
        lon_a( lon_a < -180) = lon_a(lon_a < -180)+360;
        if size(lat_a0,1) ~= size(lat_a,1); lat_a0 = lat_a0'; lon_a0 = lon_a0'; end;
        
        %check accuracy of calculation
        ind_bad = find( abs( lat_a - lat_a0) > 0.001  | ...
                        abs( lon_a - lon_a0) > 0.001 );
        if numel(ind_bad) > 0;
             disp(['Problem with Oblique Cylindrical ' ...
                    'Projection Equations!: Diff = ' ...
                    'Lat: ' num2str(max(lat_a-lat_a0)) ' Lon: ' ...
                    num2str(max(lon_a-lon_a0)) ]);
                line(ind_bad) = nan;
                sample(ind_bad) = nan;
                %return;
        end;
        
        %finally calculate line and sample
        line = line_proj_offset + lon_a.*map_res + 1;
        sample = sample_proj_offset + lat_a.*map_res + 1;
        
%         for i=1:length(lats);
%             %get the lat/lon
%             lat = lats(i);
%             lon = lons(i);
% 
%             %define XB
%             XB = [R*cos(lat)*cos(lon) R*cos(lat)*sin(lon) R*sin(lat)]';
% 
%             %calculate XA
%             XA = M*XB;
% 
%             %calculate lat_a and lon_a
%             lon_a0 = atan2( XA(2), XA(1) ).*180/pi;
%             lat_a0 = atan2( XA(3),sqrt( XA(1).^2 + XA(2).^2 ) ).*180/pi;
% 
%             %calculate alternative values for lon_a and lat_a
%             lat_a = asin( sin(lat_p)*sin(lat) + ...
%                 cos(lat_p)*cos(lat)*cos(lon-lon_p) )*180/pi;
%             val = atan2(cos(lat)*sin(lon-lon_p), ...
%                 sin(lat_p)*cos(lat)*cos(lon-lon_p) - cos(lat_p)*sin(lat));
%             lon_a = val*180/pi - rot_p*180/pi;
%             if lon_a < -180; lon_a = lon_a + 360; end;
% 
%             %check accuracy of calculation
%             if ( (abs( lat_a - lat_a0) > 0.001) | ...
%                     (abs( lon_a - lon_a0) > 0.001) );
%                 disp(['Problem with Oblique Cylindrical ' ...
%                     'Projection Equations!: Diff = ' ...
%                     'Lat: ' num2str(lat_a-lat_a0) ' Lon: ' ...
%                     num2str(lon_a-lon_a0) ]);
%                 line = nan;
%                 sample = nan;
%                 return;
%             end;
% 
%             %finally calculate line and sample
%             line(i) = line_proj_offset + lon_a.*map_res + 1;
%             sample(i) = sample_proj_offset + lat_a.*map_res + 1;
%         end;     
    case{'polar stereographic'}
%         lat_cen =pi./2;
%         lon_cen = center_lon;
%         f = 1./(1+ cos(center_lat) );
%         %lat_cen = center_lat;
%         k = f.*2./(1+sin(lat_cen).*sin(lats)+cos(lat_cen).*cos(lats).*cos(lons-lon_cen));
%         lat_cen = pi./2;
%         x = k.*cos(lats).*sin(lons-lon_cen);
%         y = k.*(cos(lat_cen).*sin(lats) - sin(lat_cen).*cos(lats).*cos(lons-lon_cen));
%         x = x.*180./pi;
%         y = y.*180./pi;
%         sample = +sample_proj_offset + soff + x.*map_res;%.*cos(center_lat);
%         line = line_proj_offset + loff - y.*map_res;%.*cos(center_lat);
        mstruct = defaultm('stereo');
        mstruct.origin = [center_lat.*180./pi ...
        center_lon.*180./pi];
        mstruct.geoid = [R 0];
        mstruct = defaultm(mstruct);
        %                    x = (samples-(2.0.*soff) - sample_proj_offset).*map_scale.*1000;
        %                    y = (1).*(line_proj_offset - lines + (1.0.*loff)).*map_scale.*1000;
        [x,y] = mfwdtran(mstruct,lats.*180./pi,lons.*180./pi);
        sample = x./map_scale./1000+1.0+sample_proj_offset;
        line = -1.*(y./map_scale./1000-1.0-line_proj_offset);         
    case{'equirectangular'};
        %convert the lat and lon polygon vectors into line and sample
        lon = lons.*180./pi;
        lat = lats.*180./pi;
        if isfield(infom,'longitudesystem');
            lonsys = (infom.longitudesystem);
            if ischar(lonsys); lonsys = str2num(lonsys); end;
            if lonsys == 180;
                lon(lon > 180) = lon(lon > 180) - 360;
            elseif lonsys == 360;
                lon(lon < 0) = lon(lon < 0) + 360;
            end;
        end;
        lat_cen = center_lat.*180./pi;
        lon_cen = center_lon.*180./pi;
        line = line_proj_offset + loff - lat.*map_res;
        sample = sample_proj_offset + soff + ...
            (lon-lon_cen).*cos(lat_cen.*pi/180).*map_res;
    case{'equirectangular cylindrical'};
        %convert the lat and lon polygon vectors into line and sample
        lon = lons.*180./pi;
        lat = lats.*180./pi;
        if isfield(infom,'longitudesystem');
            lonsys = (infom.longitudesystem);
            if ischar(lonsys); lonsys = str2num(lonsys); end;
            if lonsys == 180;
                lon(lon > 180) = lon(lon > 180) - 360;
            elseif lonsys == 360;
                lon(lon < 0) = lon(lon < 0) + 360;
            end;
        end;
        lat_cen = center_lat.*180./pi;
        lon_cen = center_lon.*180./pi;
        line = line_proj_offset + loff - lat.*map_res;
        sample = sample_proj_offset + soff + ...
            (lon-lon_cen).*cos(lat_cen.*pi/180).*map_res;
    otherwise
        disp(['Unsupported Map Projection: ' proj_type ' !!']);
        line = nan;
        sample = nan;
        return
end;
