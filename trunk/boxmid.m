function boxmid=boxmid(boxes,indo)
% boxmid=boxmid(boxes)
% boxmid=boxmid(boxes,1)
% boxmid=boxmid(boxes,2)
%
% INPUT:
% 
% boxes     a lrtb matrix suitable for FILLBOX
% indo      the index of the output coordinate if you want only one
%
% OUTPUT:
%
% boxmid    good text location with 'HorizontalAlignment','Center'
%
% Last modified by fjsimons-at-alum.mit.edu, 08/08/2014

defval('indo',0)
boxmid=[boxes(:,1)+(boxes(:,2)-boxes(:,1))/2 ...
	boxes(:,4)+(boxes(:,3)-boxes(:,4))/2];

if indo~=0
  boxmid=boxmid(indo);
end
