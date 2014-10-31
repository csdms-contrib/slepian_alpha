function kids=getkids(aks,varargin)
% kids=GETKIDS(ah)
% kids=GETKIDS(ah,ind)
%
% Get axis handle's children handles.
%
% Last modified by fjsimons-at-alum.mit.edu, June 10th, 2004

kids=get(aks,'Children');

if nargin>1
  kids=kids(varargin{1});
end
