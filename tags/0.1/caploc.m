function [lon2,lat2]=caploc(lonlat,TH,N,method)
% [lon2,lat2]=CAPLOC([lon1 lat1],TH,N,method)
%
% Calculates the location of a circle on the sphere
%
% INPUT:
%
% [lon1 lat1]      Longitude, latitude of center(s) [degrees]
% TH               Radius, in degrees [default: 15]
% N                Number of points [default: 100]
% method           1 Cartesian i.e. lon/lat [default]
%                  2 Mollweide projection
%
% OUPUT:
%
% [lon2,lat2]      The requested points [degrees or Mollweide]
%
% See also: LATITUDE, LONGITUDE
%
% Last modified by fjsimons-at-alum.mit.edu, 11/05/2010

defval('lonlat',[55 20])
defval('TH',15)

defval('N',100)
defval('method',1)
angl=linspace(0,360,N);

delta=TH*fralmanac('DegDis','Earth')/1000;

[lon2,lat2]=grcazim(lonlat,delta,angl);

if method==2
   [lon2,lat2]=mollweide(lon2*pi/180,lat2*pi/180);
end
