function [gcdkm,delta]=grcdist(lon1lat1,lon2lat2)
% [gcdkm,delta]=GRCDIST([lon1 lat1],[lon2 lat2])
%
% Calculates the distance between points on a great circle.
%
% INPUT:
%
% [lon1 lat1]  longitude and latitude of the start points [degrees]
% [lon2 lat2]  longitude and latitude of the end points [degrees]
%
% OUTPUT:
%
% gcdkm        great-circle distance [km]
% delta        great-circle distance [degrees]
%
% No ellipticity correction.
%
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2007

% Use this in conjunction with setvar.pl!

% Conversion to radians
lon1lat1=lon1lat1*pi/180;
lon2lat2=lon2lat2*pi/180;

[lon1,lat1]=deal(lon1lat1(:,1),lon1lat1(:,2));
[lon2,lat2]=deal(lon2lat2(:,1),lon2lat2(:,2));

% If this returns a complex just take the real part
dist=2*asin(sqrt((sin((lat1-lat2)/2)).^2 + ...
    cos(lat1).*cos(lat2).*(sin((lon1-lon2)/2)).^2));

gcdkm=dist*fralmanac('Radius')/1000;
delta=dist*180/pi;



