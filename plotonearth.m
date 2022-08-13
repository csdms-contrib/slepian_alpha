function varargout=plotonearth(data,conts,lon,lat)
% pc=PLOTONEARTH(data,conts,lon,lat)
%
% Maps data onto a sphere, with optional continents, and in an
% absolute longitude-latitude coordinate frame if so requested.
%
% INPUT:
%
% data      The two-dimensional data matrix to be rendered
% conts     1 continents will be plotted
%           0 continents will not be plotted
% lon       Optional longitudes for the data matrix columns [degrees]
% lat       Optional latitudes for the data matrix rows [degrees]
%
% OUTPUT:
%
% pc         The handle to the continents plotted
%
% See also PLOTONSPHERE, PLOTPLM
%
% Last modified by fjsimons-at-alum.mit.edu, 09/18/2017

% Default inputs
defval('data',rand(100,200))
defval('conts',1)
defval('lon',[])
defval('lat',[])

% Check input sizing in case you supply your own lon/lat grid
if ~isempty(lon) || ~isempty(lat)
  if size(data)~=size(lon) | size(data)~=size(lat)
    error('All input arrays need to be of equal size')
  end
end

if isempty(lon)
  % Make a relative mapping grid for the data
  [ny,nx]=deal(100);
  [lon,lat]=meshgrid(linspace(0,360,nx),linspace(90,-90,ny));
end

% Should you want the continents
if conts==1
  % Plot the continents in this three-dimensional rendering
  [~,pc]=plotcont([0 90],[360 -90],3);
end

% Make the mapping sphere  
[x,y,z]=sph2cart(lon*pi/180,lat*pi/180,ones(size(lat)));

% Render the data onto the mapping sphere
surface(x,y,z,'FaceColor','texture','Cdata',data);

% Cosmetics
axis image
shading flat
view(350,45)

% Optional outputs
varns={pc};
varargout=varns(1:nargout);
