function [dem,dax,ziro]=sergeicol
% [dem,dax,ziro]=SERGEICOL
%
% Returns a topography-type colormap and some of its diagnostics
%
% OUTPUT:
%
% dem        The color map
% dax        The data values of the extremes of the color map
% ziro       The index into the colormap where the data is zero
%
% EXAMPLE:
%
% If you plot data, and use this colorbar, you need to set the axis
% limits to the output dax or to something with the same ratio to get the
% zero crossing at the same point. If you want to change the axis limits
% to something with a different ratio, you need to interpolate again!
%
% See CAX2DEM, DEMMAP
%
% Last modified by fjsimons-at-alum.mit.edu, 09/12/2012

if ~exist('sergeim')
  try
    load(fullfile(getenv('IFILES'),'COLORMAPS','sergeim'))
  catch
    load('sergeim')
  end
end

% Resample at equal intervals - high, so discontinuity is sampled well
% Make sure to include the zero % Verify using RGBPLOT
reso=128;
demreso=nan(reso+1,3);
for index=1:3
  [demreso(:,index),req]=...
      discinter(sergeil,sergeim(:,index),...
		unique([0 linspace(min(sergeil),max(sergeil),reso)]));
end

% This is the new resampled colormap
dem=demreso;
% These are the data extremities of this colormap
dax=minmax(sergeil);
% And this is the index where the data is zero
ziro=find(req==0);

