function ext=lrtb2ext(lrtb)
% ext=LRTB2EXT(obj)
%
% Transforms (left-right-top-bottom) 'edge positions' of a graphics
% object into an equivalent box 'extent' (left-bottom-width-height)

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
% Last modified by fjsimons-at-alum.mit.edu, 05/26/2021

ext=[lrtb(1) lrtb(4) lrtb(2)-lrtb(1) abs(lrtb(4)-lrtb(3))];
