function g=grey(val)
% g=grey(val)
%
% Will return RGB for grey values (0-10).
%
% INPUT:
%
% val      0 (black) to 10 (white) [default: 7]
%
% OUTPU:
%
% g        An RGB vector with the grey scale
%
% Last modified by fjsimons-at-alum.mit.edu, 10/06/2006

defval('val',7)

if length(val)>1
  for index=1:length(val(:))
    g(index,1:3)=[1 1 1]*val(index)/10;
  end
else
  g=[1 1 1]*val/10;
end

