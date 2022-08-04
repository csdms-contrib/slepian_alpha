function varargout=randinpolyradius(dom,N,satrad,varalt,rotcoords)
% [rad,theta,phi]=randinpolyradius(dom,N,satrad,varalt,rotcoords) 
%
% Generates N random points within the region dom (like randinpoly) and 
% also within a given radius band. This is useful to generating synthetic data
% at a specific satellite altitude.
%
% INPUT:
%
% dom          Named region such as eurasia, namerica, samerica, etc. OR
%                Radius of the concentration region (degrees) OR
%                Two Radii (outer then inner) of the concentration region 
%                  which will be a donut (degrees)
% N            Number of points to generate. The resulting number will be
%                close but not exact
% satrad       The mean satellite radius for your data
% varalt       The full range of the radius interval for your data. You will
%                get data from satrad+varalt/2 to satrad-varalt/2.
% rotcoords    Euler angles in case you would like rotated data points 
%                (such as if you made a spherical cap. see rottp)
%
% OUTPUT:
%
% rad	    Radius of the random points
% theta     Co-latitudes of the random points (radians)
% phi       Longitude of the random points (radians)
%
% EXAMPLE:
%
% randinpolyradius('demo1')
% 
% Last modified by charig@arizona.edu, 6/16/2022
% Last modified by plattner-at-alumni.ethz.ch, 6/27/2016


defval('dom','namerica')
defval('N',200)
defval('satrad',400)
defval('varalt',100)
defval('rotcoords',[])


% Now generate random points within target region

if ischar(dom) && ~strncmp(dom,'demo',4)
    % a region
    [lon,lat]=randinpoly(dom,N);
    % Transform latitude into colatitude:
    cola=90-lat;

elseif length(dom)==1
    % a spherical cap
    [lon lat]=randpatch(N,max(dom),0,0);
    % Transform latitude into colatitude:
    cola=90-lat;

elseif length(dom)==2
    % a donut
    [lon lat]=randpatch(N,max(dom),0,0);
    % Transform latitude into colatitude:
    cola=90-lat;
    % Take out "inner cap"
    index=find(cola>=min(dom));
    lon=lon(index);
    lat=lat(index);
    cola=cola(index);

elseif strcmp(dom,'demo1')
    [rad,theta,phi]=randinpolyradius('namerica',1000,400,200);

    XY = namerica();
    figure
    plot(XY(:,1),XY(:,2),'k-'); axis equal; grid on
    hold on
    scatter(phi*180/pi,90-theta*180/pi,10,rad,'filled')
    return

end             
        
% And give them random altitudes around satrad
rad=satrad + (rand(length(lon),1)-0.5)*varalt;
fprintf('Number of random points created: %d\n',length(lon))
% Get theta and phi, which is in radians:
theta=cola*pi/180;
phi=lon*pi/180;

% And rotate the points to the right location
if ~isempty(rotcoords)
    [theta,phi]=rottp(theta,phi,0,-rotcoords(2)*pi/180,-rotcoords(1)*pi/180);
end
    
   
% Prepare output
vars={rad,theta,phi};
varargout=vars(1:nargout);
    
    
