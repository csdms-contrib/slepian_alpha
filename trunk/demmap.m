function varargout=demmap
% DEMMAP
% [dem,dax,ziro]=DEMMAP;
% 
% Sets the colormap to one good for digital elevation models, or merely
% outputs three variables useful to do this elsewhere.
%
% OUTPUT:
%
% dem        The color map
% dax        The data values of the extremes of the color map
% ziro       The index into the colormap where the data is zero
%
% See also CAX2DEM, SERGEICOL, ADDCB
%
% Last modified by fjsimons-at-alum.mit.edu, 09/29/2008

if ~exist('demmapm')
  load(fullfile(getenv('IFILES'),'COLORMAPS','topo'))
end

% Supply the minimum and maximum data values to which these apply
% load topo; mima=minmax(topo);
mima=[-7473 5731];

if nargout
  varargout{1}=demmapm;
  varargout{2}=mima;
  varargout{3}=36; % Used to be 37 for some reason
else
  colormap(demmapm)
  caxis(mima)
end

