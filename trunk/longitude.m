function [x,y,z]=longitude(lon,meth,nums)
% [x,y,z]=LONGITUDE(lon,meth,nums)
%
% Yields the coordinates of a line of constant longitude.
%
% INPUT:
%
% lat       Constant longitude, in degrees
% meth      1 3-D coordinates
%           2 2-D coordinates
%           3 2-D Mollweide coordinates
% nums      A scalar number with the number of points
%
% OUPUT:
%
% x,y,z     Coordinates fit to plot outside of the unit sphere
%
% See also: LATITUDE, CAPLOC
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 16th, 2004

defval('lat',20)
defval('meth',1)
defval('nums',100)

lat=linspace(-90,90,nums);

switch meth
 case 1
  [x,y,z]=sph2cart(lon/180*pi,lat/180*pi,1.01);
 case 2
  x=repmat(lon,size(lat));
  y=lat;
  z=[];
 case 3
  warning off
  [x,y]=mollweide(lon*pi/180,lat*pi/180,pi);
  z=[];
  warning on
end
