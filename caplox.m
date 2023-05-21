function [lon2,lat2]=caplox(lonlat,TH,N,method)
% [lon2,lat2]=CAPLOX([lon1 lat1],TH,N)
%
% Calculates the location of a circle on the sphere crossed with the
% outlines of the continents.
%
% INPUT:
%
% [lon1 lat1]      Longitude, latitude of center(s) [degrees]
% TH               Radius, in degrees [default: 15]
% N                Number of points [default: 100]
%
% OUPUT:
%
% [lon2,lat2]      The requested points [degrees or Mollweide]
%
% See also: LATITUDE, LONGITUDE, CAPLOC
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2021

defval('lonlat',[-25 10])
defval('TH',15)

defval('N',100)
defval('method',1)
angl=linspace(0,360,N);

delta=TH*fralmanac('DegDis','Earth')/1000;

[lon2,lat2]=grcazim(lonlat,delta,angl);

figure(gcf)
clf
plotcont([],[],11,[],[],lonlat)
hold on


keyboard
hold off

