function nolabels(axh,w)
% NOLABELS(handels) % for x and y 
% NOLABELS(handels,1) % For x
% NOLABELS(handels,2) % for y
% NOLABELS(handels,3) % for x and y
%
% SEE ALSO: 
%
% NOTICKS
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/08/2015

defval('axh',gca)
defval('w',3)

switch w 
  case 3
    set(axh,'XTickLabel',[],'YTickLabel',[])
  case 2
    set(axh,'YTickLabel',[])
  case 1
    set(axh,'XTickLabel',[])
end

