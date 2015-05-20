function boxmid=boxmid(boxes,indo)
% boxmid=boxmid(boxes)
% boxmid=boxmid(boxes,1)
% boxmid=boxmid(boxes,2)
%
% 'boxes' being a lrtb matrix good for fillbox
% 'boxmid' is a good text location with 'HorizontalAlignment','Center'
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 11th, 2003

defval('indo',0)
boxmid=[boxes(:,1)+(boxes(:,2)-boxes(:,1))/2 ...
	boxes(:,4)+(boxes(:,3)-boxes(:,4))/2];

if indo~=0
  boxmid=boxmid(indo);
end
