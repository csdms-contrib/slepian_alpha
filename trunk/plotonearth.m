function varargout=plotonearth(data,conts,lon,lat)
% PLOTONEARTH(data)
% PLOTONEARTH(data,1) with continents
% PLOTONEARTH(data,1,lon,lat) if it matters where
%
% h=PLOTONEARTH returns axis handle to the continents
%
% Plots data on a globe, with continents or not.
% Without relief but with optional continents, and with
% the possibility of an absolute coordinate frame. 
%
% See also PLOTONSPHERE, PLOTPLM
%
% Last modified by fjsimons-at-alum.mit.edu, June 26th, 2003

defval('conts',1)
defval('lon',[])
defval('lat',[])

if ~isempty(lon)
  if size(data)~=size(lon) | size(data)~=size(lat)
    error('Wrong size arrays')
  end
end

if conts==1
  % Plot the continents
  [jk1,jk2,cont]=plotcont([0 90],[360 -90]);
  delete(jk2)
  lonc=cont(:,1)/180*pi;
  latc=cont(:,2)/180*pi;
  rad=repmat(1.001,size(latc));
  [xx,yy,zz]=sph2cart(lonc,latc,rad);
  pc=plot3(xx,yy,zz,'k-','LineWidth',1.5);
  hold on
end

if isempty(lon)
  % Make sphere for the data
  lons=linspace(0,360,100)/180*pi;
  lats=linspace(90,-90,100)/180*pi;
  [lons,lats]=meshgrid(lons,lats);
else
  lons=lon;
  lats=lat;
end

rads=ones(size(lats));

[x,y,z]=sph2cart(lons,lats,rads);

surface(x,y,z,'FaceColor','texture','Cdata',data);

axis image
shading flat
view(140,30)

nargs={'pc'};
for ind=1:nargout
  varargout{ind}=eval(nargs{ind});
end
