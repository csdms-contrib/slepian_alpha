function [flon,tlat]=project(lon,lat,c11,cmn)
% [flon,tlat]=project(lon,lat,c11,cmn)
%
% Give it true latitudes and longitudes and projects coordinates
% just like EQUISTAT does.
%
% See also UNPROJECT, EQUISTAT, POLARGRID
%
% EXAMPLE I 
%
% C11=[112.6989 -5.7875];
% CMN=[156.9554 -48.2125];
% [on,tw,XY]=plotcont([90 10],[180 -60]); hold on
% [flon,tlat]=project(XY(:,1),XY(:,2),C11,CMN);
% plot(flon,tlat,'r')
% [lon,lat]=unproject(flon,tlat,C11,CMN);
% plot(lon,lat,'g')
%
% EXAMPLE II
%
% C11=[112.6989 -5.7875];
% CMN=[156.9554 -48.2125];
% ddir= '/home/fjsimons/MyPapers/2003/JGR-2003/DATA/';
% load(fullfile(ddir,'tdint'))
% imagef(C11,CMN,tdint); axis image
% cb=cax2dem([-8000 1500]);
% [on,tw,XY]=plotcont([90 10],[180 -60]); delete(tw)
% hold on
% [flon,tlat]=project(XY(:,1),XY(:,2),C11,CMN);
% plot(flon,tlat,'r')
%
% Last modified by fjsimons-at-alum.mit.edu, 09/24/2008

XIM=[c11(1) cmn(1)];
YIM=[c11(2) cmn(2)];
 
tlat=lat;
flon=c11(1)+(lon-c11(1)).*cos(lat*pi/180);
shift=[cmn(1)-c11(1)]*(1-cos(tlat*pi/180))/2;
flon=flon+shift;


