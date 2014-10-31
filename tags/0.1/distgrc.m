function [lon,lat]=distgrc(lon1lat1,lon2lat2,dist)     
% [lon,lat]=distgrc([lon1 lat1],[lon2 lat2],dist)
%
% Calculates coordinates of points at a certain distance on a great circle
%
% INPUT:
%
% [lon1 lat1]  longitude and latitude of the start point [radians]
% [lon2 lat2]  longitude and latitude of the end point [radians]
% dist         distance of the requested points [radians]
%
% OUTPUT:
%
% [lon lat]    longitude and latitude of the requested points [radians]
%
% EXAMPLE:
%
% lon1lat1=[110 12]; lon2lat2=[220 -23]; fax=1000/fralmanac('Radius');
% dist=linspace(0,grcdist(lon1lat1,lon2lat2)*fax,20);
% [lon,lat]=distgrc(lon1lat1*pi/180,lon2lat2*pi/180,dist);
% plot(lon1lat1(:,1),lon1lat1(:,2),'s'); hold on
% plot(lon2lat2(:,1),lon2lat2(:,2),'s')
% plot((lon+(lon<0)*2*pi)*180/pi,lat*180/pi,'+'); hold off
%
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2007

[lon1,lat1]=deal(lon1lat1(1),lon1lat1(2));
[lon2,lat2]=deal(lon2lat2(1),lon2lat2(2));

% Calculates initial course
tc=truecourse([lon1 lat1]*180/pi,[lon2 lat2]*180/pi)*pi/180;

lat =asin(sin(lat1)*cos(dist)+cos(lat1)*sin(dist)*cos(tc));
dlon=atan2(sin(tc)*sin(dist)*cos(lat1),cos(dist)-sin(lat1)*sin(lat));
lon=mod(lon1-dlon +pi,2*pi )-pi;

% Output in column vectors
lon=lon(:);
lat=lat(:);
