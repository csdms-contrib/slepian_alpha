function varargout=fridplot(lonlon,latlat,varargin)
% ph=FRIDPLOT(lonlon,latlat,'Property','Value')
%
% Plots a Cartesian grid
%
% INPUT:
%
% lonlon       Longitudes, as a vector or from MESHGRID
% latlat       Latitudes, as a vector or from MESHGRID
% 'Property'   A list of handle properties
% 'Value   '   A list of handle property values
%
% OUTPUT:
%
% ph           Two handles to the plotted graphics objects
% 
% EXAMPLE:
%
% [lonlon,latlat]=meshgrid(linspace(10,20,25),linspace(40,60,15));
% fridplot(lonlon,latlat); openup(gca,5); openup(gca,6)
%
% [lonlon,latlat]=equistat([10 60],[20 40],25,15);
% fridplot(lonlon,latlat); axis([-2 10 -22 2])
%
% Last modified by fjsimons-at-mit.edu, 08/08/2008

if min(size(lonlon))==1 && min(size(latlat))==1
  [lonlon,latlat]=meshgrid(lonlon,latlat);
end

% Produce lift-the-pen points
latlat=adrc(latlat,NaN);
lonlon=adrc(lonlon,NaN);

% Do the plotting
ph=plot(lonlon(:),latlat(:),'k',indeks(lonlon',':'),indeks(latlat',':'),'k');

% Cosmetic adjustments
if nargin>2; set(ph,varargin{1:end}); end

% Prepare output
varn={ph}; varargout=varn(1:nargout);


