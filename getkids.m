function kids=getkids(ah,varargin)
% kids=GETKIDS(ah,ind)
%
% Get axis handle's children handles.
%
% INPUT:
%
% ah       An axis handle
% ind      An index of you want a particular one
%
% Last modified by fjsimons-at-alum.mit.edu, 05/26/2021

kids=get(ah,'Children');

% Subselection if you know what you want
if nargin>1
  kids=kids(varargin{1});
end
