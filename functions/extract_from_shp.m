function [x_center,y_center,x_leftedge,y_leftedge,x_rightedge,y_rightedge] = extract_from_shp(feature_name,planet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract values for mapped banks and centerline from shapefiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   feature_name : unique river name identifier made during mapping (e.g. BigBend)
%   planet       : name of planet (e.g., 'Titan')
%   
% OUTPUT
%   x_center     : longitude (deg) of mapped centerline points
%   y_center     : latitude (deg) of mapped centerline points
%   x_leftedge   : longitude (deg) of mapped left bank points
%   y_leftedge   : latitude (deg) of mapped left bank points
%   x_rightedge  : longitude (deg) of mapped right bank points
%   y_rightedge  : latitude (deg) of mapped right bank points
%
% Note, the name right bank and left bank is arbitrary and was decided during 
% mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    if strcmp(planet,'Titan')
        % load centerline
        feature_shp = feature_name+"CenterLine.shp"; 
        % Load feature edge lines
        leftbank_shp = feature_name+"LeftEdge.shp"; 
        rightbank_shp = feature_name+"RightEdge.shp"; 
    elseif strcmp(planet,'Earth')
        feature_shp = fullfile(feature_name,'centerline.shp');
        leftbank_shp = fullfile(feature_name,'leftEdge.shp');
        rightbank_shp = fullfile(feature_name,'rightEdge.shp');
    else
        error('planet name required in function call')
    end
    
    if ~exist(feature_shp,'file')
        %fprintf('no shape file available for %s\n',feature_shp)
        x_center = [];
        y_center = [];
        x_leftedge = [];
        y_leftedge = [];
        x_rightedge = [];
        y_rightedge = [];
        return
    end
    
    warning('off','map:shapefile:buildingIndexFromSHP'); 
    fs = shaperead(feature_shp);
    
    % Check for and get rid of nan nodes
    x_center = fs.X(~isnan(fs.X))';
    y_center = fs.Y(~isnan(fs.Y))';
    
    leftbank = shaperead(leftbank_shp);
    rightbank = shaperead(rightbank_shp);
    
    if strcmp(feature_name,'EquatorialEx')
        % Note from Frank:
        % Shapefiles of EquatorialEx have a mistake, and are subdivided into
        % separate lines. Use this code to load in instead.
        % x1 = [edge1(1).X,edge1(2).X,edge1(3).X,edge1(4).X];
        % y1 = [edge1(1).Y,edge1(2).Y,edge1(3).Y,edge1(4).Y];
        % x2 = [edge2(1).X,edge2(2).X];
        % y2 = [edge2(1).Y,edge2(2).Y];
        % x_leftedge = [leftbank(1).X,leftbank(2).X,leftbank(3).X,leftbank(4).X];
        %  y_leftedge = [leftbank(1).Y,leftbank(2).Y,leftbank(3).Y,leftbank(4).Y];
        %  x_rightedge = [rightbank(1).X,rightbank(2).X,rightbank(3).X,rightbank(4).X];
        %  y_rightedge = [rightbank(1).Y,rightbank(2).Y,rightbank(3).Y,rightbank(4).Y];
        
        x_leftedge = leftbank.X;
        y_leftedge = leftbank.Y;
        x_rightedge = rightbank.X;
        y_rightedge = rightbank.Y;
        
        x_leftedge = rmmissing(x_leftedge);
        y_leftedge = rmmissing(y_leftedge);
        x_rightedge = rmmissing(x_rightedge);
        y_rightedge = rmmissing(y_rightedge);
    else
        x_leftedge = rmmissing(leftbank.X);
        y_leftedge = rmmissing(leftbank.Y);
        x_rightedge = rmmissing(rightbank.X);
        y_rightedge = rmmissing(rightbank.Y);
    end
    
    
    x_leftedge = x_leftedge(:);
    y_leftedge = y_leftedge(:);
    x_rightedge = x_rightedge(:);
    y_rightedge = y_rightedge(:);

end