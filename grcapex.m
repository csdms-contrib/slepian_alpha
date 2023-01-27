function [latmax,lonmax]=grcapex(lat1,lon1,lat2,lon2)
% [latmax,lonmax]=grcapex(lat1,lon1,lat2,lon2)
%
% Calculates highest points on a great circle
%
% INPUT:
%
% lat1,lon1       A geographical point, given in radians
% lat2,lon2       Another geographical point, in radians
%
% OUTPUT:
%
% latmax,lonmax   The position of the highest point on the great
%                 circle joining the two input positions, in radians
%
% Last modified by fjsimons-at-alum.mit.edu, 01/26/2023

% Calculates the great circle distance between two points, in radians
% d = acos(sin(lat1).*sin(lat2)+cos(lat1).*cos(lat2).*cos(lon1-lon2));
d=grcdist([lon1 lat1]*180/pi,[lon2 lat2]*180/pi)/6371;

% Calculates the true course at the point d on the line
tc1 = acos((sin(lat2)-sin(lat1).*cos(d))./(sin(d).*cos(lat1)));

%  Highest latitude reached
latmax=acos(abs(sin(tc1).*cos(lat1)));

% Calculate longitude given latitude on great circle (input in degrees)
l12=lon1-lon2;
A =sin(lat1).*cos(lat2).*cos(latmax).*sin(l12);
B=sin(lat1).*cos(lat2).*cos(latmax).*cos(l12)-cos(lat1).*sin(lat2).*cos(latmax);
lonex=atan2(B,A);
lonmax=lon1+lonex+(sin(l12)<0)*pi;
lonmax=rem(lonmax,2*pi)+(lonmax<0)*2*pi;
