function [cb,cm]=cax2dem(cax,ori,lronly)
% [cb,cm]=CAX2DEM(cax,ori,lronly)
%
% Maps the limits of the data into a nice discontinuous colormap.
% This just resets the colormap; no direct indexing of data
% is necessary as opposed to JOINCOLMAP.
%
% INPUT:
%
% cax       The low/high data at which the colormap saturates
% ori       'hor' or 'ver' [default] position of the colorbar
% lronly    If you're only needing the left or right color portion
%
% OUTPUT:
%
% cb        The handle to the color bar
% cm        The actual color map
%
% EXAMPLE: for topography of Australia:
%
% cb=cax2dem([-7000 1500]);
%
% EXAMPLE:
%
% frsdemo1('Earth'); frsdemo1('Venus'); frsdemo1('Moon'') etc
%
% See also: ADDCB, JOINCOLMAP, SERGEICOL, DEMMAP
%
% Last modified by fjsimons-at-alum.mit.edu, 01/21/2022

defval('ori','ver')
defval('lronly','all')

% Get the color map and its attributes
[dem,dax,ziro]=sergeicol;

% So it "knows" the data values belonging to the original colors
if cax(2)>dax(2) || cax(1)<dax(1)
  warning([ 'The requested data axis range must be '...
      'smaller than the limits of the colormap.'])
end

% The zero value in this color map comes at 'ziro'; 
% this is the same ratio between high and low in the real data
ratd=abs(dax(2)/dax(1));

% What is the ratio for our desired caxis value?
rat=abs(cax(2)/cax(1));

% Length of the current color map
len=size(dem,1);

% Add colors to the color map - on either side, should be fine; all
% you're doing is making sure the number of colors is right and Matlab
% interpolates the rest anyway
switch lronly
 case 'left'
  cm=dem(1:ziro,:);
 case 'right'
  cm=dem(ziro+1:len,:);
 otherwise
  cm=interp1(1:len,dem,...
	     [linspace(1,ziro,round((len-ziro)/rat)) ziro+1:len]);
end

% Apply the choices
colormap(cm)
caxis(cax)

% Default to vertical orientation
cb=colorbar(ori);
