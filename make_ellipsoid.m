function ellipsoid = make_ellipsoid(planetname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a simple ellipsoid for the planet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   planetname : name of planet e.g., 'Titan'
% OUTPUT
%   ellipsoid  : simple ellipsoid to describe planet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if strcmp(planetname,'Titan')
        planet_r = 2575000;   % avg radius meters
    end
    
    if strcmp(planetname,'Earth')
        planet_r = 6378000;   % avg radius meters
    end
    
    ellipsoid = referenceEllipsoid;
    ellipsoid.Name = planetname;
    ellipsoid.LengthUnit = 'meter';
    ellipsoid.SemimajorAxis = planet_r;
    ellipsoid.SemiminorAxis = planet_r; % make it a sphere


end