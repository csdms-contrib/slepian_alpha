function noticks(varargin)
% NOTICKS(handels)
% NOTICKS(handels,1)
% NOTICKS(handels,2)
% NOTICKS(handels,3)
%
% Last modified by fjsimons-at-alum.mit.edu, 11/13/2013

defval('handels',gca)

poshan={'Xtick','YTick','ZTick'};

if nargin==1
  set(varargin{1},'XTick',[],'YTick',[],'ZTick',[])
else
  set(varargin{1},poshan{varargin{2}},[])
end
