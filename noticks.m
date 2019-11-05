function noticks(varargin)
% NOTICKS(handels) % for x, y and z
% NOTICKS(handels,1) % for x
% NOTICKS(handels,2) % for y
% NOTICKS(handels,3) % for z
%
% SEE ALSO: 
%
% NOLABELS
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/08/2015

defval('handels',gca)

poshan={'Xtick','YTick','ZTick'};

if nargin==1
  set(varargin{1},'XTick',[],'YTick',[],'ZTick',[])
else
  set(varargin{1},poshan{varargin{2}},[])
end

