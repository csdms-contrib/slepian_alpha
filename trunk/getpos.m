function pos=getpos(varargin)
% pos=GETPOS(ah)
% pos=GETPOS(ah,ind)
%
% Gets axis handle position or just one of its four elements. 
% If input is 1xN vector returns Nx4 or Nx3 positions.
%
% INPUT:
%
% ah      The axis handle whose position you want
% ind     The index in the position vector you want
%
% OUTPUT:
%
% pos     The position vector
%
% See also GETTIT
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

defval('ah',gca)

pos=get(varargin{1},'Position');

if iscell(pos)
  vr=length(pos{1});
  pos=reshape([pos{:}],vr',length([pos{:}])/vr')';
end

if nargin>1
  pos=pos(:,varargin{2});
end

