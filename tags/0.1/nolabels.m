function nolabels(axh,w)
% NOLABELS(handels)
% NOLABELS(handels,1)
% NOLABELS(handels,2)
%
% Last modified by fjsimons-at-mit.edu, Nov. 2nd, 2001

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

