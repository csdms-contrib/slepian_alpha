function [lon,lat]=randinpoly(dom,N)
% [lon,lat]=randinpoly(dom,N) 
%
% Generates less than N equal-area random points within the region dom 
%
% INPUT:
%
% dom 		Named region such as eurasia, namerica, samerica, etc.
% N 		Number of points to generate. The resulting number will be less!
%
% OUTPUT:
%
% lon		longitudes of the random points
% lat 		latitudes of the random points
% 
% Last modified by plattner-at-alumni.ethz.ch, 6/27/2016

XY=eval(sprintf('%s()',dom));
% Use randpatch to make circle with equal area random point

center=[sum(XY(:,1)) sum(XY(:,2))]/length(XY(:,1));
distx=XY(:,1)-center(1);
disty=XY(:,2)-center(2);
allr=sqrt(distx.^2+disty.^2);
r=max(allr);
[lon,lat]=randpatch(N,r,center(1),90-center(2));
% Sometimes the XY long goes beyond 360. Make lon go beyond by the same
% amount
if max(XY(:,1))>360
  shift=max(XY(:,1))-mod(max(XY(:,1)),360);
  lon2=lon+shift;
  lon=[lon;lon2];
  lat=[lat;lat];
end

% Now delete all points outside the region
index=inpolygon(lon,lat,XY(:,1),XY(:,2));
lon=lon(index);
lat=lat(index);

fprintf('More or less randomly generated %d points in %s\n',...
	sum(index),dom);

