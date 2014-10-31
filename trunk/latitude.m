function [x,y,z]=latitude(lat,meth,nums)
% [x,y,z]=LATITUDE(lat,meth,nums)
%
% Yields the coordinates of a line of constant latitude.
%
% INPUT:
%
% lat       Constant latitude, in degrees
% meth      1 3-D coordinates on a sphere of radius 1.01 [default]
%           2 2-D coordinates
%           3 2-D Mollweide coordinates
% nums      A scalar number with the number of points [default 100]
%
% See also: LONGITUDE, CAPLOC
%
% Last modified by fjsimons-at-alum.mit.edu, August 4th, 2004

defval('lat',20)
defval('meth',1)
defval('nums',100)

lon=linspace(0,360,nums);

switch meth
 case 1
  [x,y,z]=sph2cart(lon/180*pi,lat/180*pi,1.01);
  z=repmat(z,size(x));
 case 2
  x=lon;
  y=repmat(lat,size(x));
  z=[];
 case 3
  warning off
  [x,y]=mollweide(lon*pi/180,lat*pi/180,pi);
  y=repmat(y,size(x));
  z=[];
  warning on
end






