function [lon,lat]=unproject(flon,tlat,c11,cmn)
% [lon,lat]=UNPROJECT(flon,tlat,c11,cmn)
%
% Give it true latitudes and fake longitudes (i.e. the x coordinates in the
% projected map) and it calculates the real [lon,lat]. With respect to the
% top latitude of the data (see EQUISTAT) and the [cmn(1)-c11(1)]/2 (see EQUISTAT).
%
% flon and flat can be just points - no meshgrid necessary
%
% See EQUISTAT, PROJECT
%
% Last modified by fjsimons-at-mit.edu, June 9th, 2001

lat=tlat;
shift=[cmn(1)-c11(1)]*(1-cos(tlat*pi/180))/2;
lon=flon-shift;
lon=c11(1)+(lon-c11(1))./cos(tlat*pi/180);




