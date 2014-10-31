function [lonlon,latlat,intlon,intlat,orlon,orlat]=equistat(c11,cmn,lonnum,latnum)
% [lonlon,latlat,intlon,intlat,orlon,orlat]=EQUISTAT(c11,cmn,lonnum,latnum)
%
% Expresses a geographical grid of regularly spaced latitudes and
% longitudes in terms of the true distances along the surface of the
% sphere. Suitable to interpolate a geographical grid to a locally
% regular Cartesian flat grid.
%
% INPUT:
%
% C11         [lon lat] of the upper left corner of the map [degrees]
% CMN         [lon lat] of the lower right corner of the map [degrees]
% lonnum      Number of samples across (E-W)
% latnum      Number of samples down (N-S)
%
% OUTPUT:
%
% lonlon      The Cartesion x-coordinates of the geographical [lon lat] grid
% latlat      The Cartesion y-coordinates of the geographical [lon lat] grid
% intlon      Regular Cartesian longitudes for interpolation using GRIDDATA
% intlat      Regular Cartesian latitudes for interpolation using GRIDDATA
%             All of the above are with respect to a top left corner [0,0]
% orlon       The longitudes of the original geographical grid
% orlat       The latitudes of the original geographical grid
%             The above two with respect to the actual top left corner 
%
% EXAMPLE I:
%
% XIM=[115 155]; YIM=[-5 -50];
% [lonlon,latlat,intlon,intlat,orlon,orlat]=...
%    equistat([XIM(1) YIM(1)],[XIM(2) YIM(2)],20,20);
% fridplot(lonlon,latlat); hold on
% co=fridplot(intlon,intlat); set(co,'Color','b') ; hold off
%
% EXAMPLE II:
%
% XIM=[100 160]; YIM=[1 -50];
%% Get the data on a regular geographical grid (longitude and latitude)
% z=flipud(etopo(fullfile(getenv('IFILES'),'TOPOGRAPHY','EARTH'),5,sort(YIM),XIM));
% subplot(121)
% imagef([XIM(1) YIM(1)],[XIM(2) YIM(2)],z); axis image; 
%% Plot the continental outlines, in geographical coordinates, on top
% [a,b,XY]=plotcont([XIM(1) YIM(1)],[XIM(2) YIM(2)]); axis tight
% cb=cax2dem([-8000 1500]); delete(cb)
%% Find interpolation points on a regular Cartesian grid (E-W and N-S)
%% contained in the original grid
% [lonlon,latlat,intlon,intlat]=...
%    equistat([XIM(1) YIM(1)],[XIM(2) YIM(2)],size(z,2),size(z,1));
%% Interpolate the geographical grid to this regular Cartesian grid
% zi=griddata(lonlon,latlat,z,intlon,intlat);
% subplot(122)
% Plot the interpolated image on the regular Cartesian grid
% imagef([0 0],[range(intlon(:)) -range(intlat(:))],zi); axis image
%% Project the continental outlines according to the same scheme
% [flon,tlat]=project(XY(:,1),XY(:,2),[XIM(1) YIM(1)],[XIM(2) YIM(2)]);
% cb=cax2dem([-8000 1500]); delete(cb)
% The next line needs to be fixed !!!!
%% hold on ; d=plot(flon...,tlat...,'k'); hold off
%
% See also UNPROJECT, PROJECT, POLARGRID, ELL2CAR, MASTERSET
%
% Last modified by fjsimons-at-alum.mit.edu, 09/24/2008

XIM=[c11(1) cmn(1)];
YIM=[c11(2) cmn(2)];

longrid=linspace(XIM(1),XIM(2),lonnum);
latgrid=linspace(YIM(1),YIM(2),latnum);

[orlon,orlat]=meshgrid(longrid,latgrid);

lonlon=(orlon-c11(1)).*cos(orlat*pi/180);
latlat=(orlat-c11(2));

% Origin in middle of data set
shift=repmat((lonlon(1,lonnum)-lonlon(:,lonnum))/2,1,lonnum);
lonlon=lonlon+shift;
latlat=latlat;

intval=max(lonlon(:,lonnum)-lonlon(:,lonnum-1));
minval=max(lonlon(:,1)); 
maxval=min(lonlon(:,lonnum)); 
intlon=repmat(minval:intval:maxval,latnum,1);
intlat=repmat(latlat(:,1),1,size(intlon,2));











