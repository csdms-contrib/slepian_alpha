function ext=lrtb2ext(lrtb)
% ext=LRTB2EXT(obj)
%
% INPUT:
%
% lrtb      Coordinates (left-right-top-bottom)
%
% OUTPUT:
% 
% ext       Extent (left-bottom-width-height)
%
% See also FILLBOX, ALLEN1, EXT2LRTB
% 
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

ext=[lrtb(1) lrtb(4) lrtb(2)-lrtb(1) abs(lrtb(4)-lrtb(3))];


