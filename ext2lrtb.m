function lrtb=ext2lrtb(obj,wf,hf)
% lrtb=EXT2LRTB(obj,wf,hf)
%
% Transforms the extent (left-bottom-width-height) of a graphics object
% into a lrtb (left-right-top-bottom) set of coordinates, e.g. to provide
% a white box under a line of text in a figure.
%
% INPUT:
%
% obj       Object handle, e.g. to a text object
% wf        Scale factor to be applied to the box width [default: 1]
%           This involves the width counting from the left.
% hf        Scale factor to be applied to the box height [default: 1]
%           This involves the height counting from the middle of the box
% INPUT:
%
% lrtb      Coordinates in the left-right-top-bottom convention
%
% See also FILLBOX, LRTB2EXT
% 
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2012

defval('wf',1)
defval('hf',1)

try
  ext=get(obj,'extent');
catch
  ext=get(obj,'position');
end

% This seems to depend on the orientation of the axis... give this a try
if strcmp(get(gca,'ydir'),'normal')
  lrtb=[ext(1) ext(1)+ext(3)*wf ext(2)+ext(4)*hf ext(2)+(ext(4)*(1-wf)/2)];
elseif  strcmp(get(gca,'ydir'),'reverse')
  lrtb=[ext(1) ext(1)+ext(3)*wf ext(2)-ext(4)*hf ext(2)+(ext(4)*(1-wf)/2)];
end
