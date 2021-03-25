function nolabels(axh,w)
% NOLABELS(axh,w)
% 
% INPUT:
%
% axh    Axis handles (default: gca)
% w      1 for x-axes
%        2 for y-axes
%        3 for x and y-axes (default)
%        4 for x, y and z-axes
%
% SEE ALSO: 
%
% NOTICKS
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/24/2021

defval('axh',gca)
defval('w',3)

switch w 
  case 3
    set(axh,'XTickLabel',[],'YTickLabel',[])
  case 2
    set(axh,'YTickLabel',[])
  case 1
    set(axh,'XTickLabel',[])
  case 4
    set(axh,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
end

