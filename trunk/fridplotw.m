function varargout=fridplotw(N,J,xyscale,varargin)
% [ph,W]=FRIDPLOTW(N,J,xyscale,'Property','Value')
%
% Plots a multiresolution wavelet decomposition grid
%
% INPUT:
%
% N            The power of the dyadic subdivision [defaulted]
% J            Maximum scale in both directions
% xyscale      Where the grid should end up for real, see PLOTONCUBE
% 
% 'Property'   A list of handle properties
% 'Value   '   A list of handle property values
%
% OUTPUT:
%
% ph           Handles to the plotted graphics objects
% W            A matrix with the coordinates of the square midpoints
% 
% EXAMPLE:
%
% N=8; J=3; f=fridplotw(N,J);
% 
% Last modified by fjsimons-at-mit.edu, 1/26/2011

% Here are all the intercepts
w=2.^(N-[0:J+1]);
rims=[0.5 sort(w(1:end-1)+0.5)];
grids=[pauli(rims,2) ; [0.5*ones(1,J+1) ; rims(2:end)]'];

defval('xyscale',[])

for in=1:size(grids,1)      
  if isempty(xyscale)
    ph{in}=fridplot(grids(in,:),grids(in,:));
  else
    ph{in}=fridplot(xyscale(1)+grids(in,:)/max(grids(:))*[xyscale(2)-xyscale(1)],...
                    xyscale(4)-grids(in,:)/max(grids(:))*[xyscale(4)-xyscale(3)]);
  end
  hold on
end

if nargout==2
  x=w(3:end);
  v=[w(2:end-1)+x w(end)];
  W=[v x v(1:end-1) ; v v(1:end-1) x];
else
  W=NaN;
end

% Cosmetic adjustments
if nargin>3
  set(cat(1,ph{:}),varargin{1:end})
end

% Prepare output
varn={ph,W}; 
varargout=varn(1:nargout);
