function varargout=supertit(ah,tstr,fs)
% [H,ta]=SUPERTIT(ah,tstr,fs)
%
% Takes a set of axis handles and and puts a title string spanning over the
% panels in question. These axis handles need to refer to panels on the same
% window row of a figure, and they need to have the same height.
%
% INPUT:
%
% ah        Axis handles
% tstr      Title string
% fs        Fontsize [default: 15]
%
% OUTPUT:
%
% H         Handle to the supertitle created
% ta        Handle to the invisible axis containing it
%
% EXAMPLE:
%
% a=krijetem([231 232 233]);
% b=supertit(a([1 2]),'This is the supertitle');
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 09/19/2016

defval('fs',15);

pos=get(ah,'position');
if iscell(pos)
  pos=cat(1,pos{:});
end
left=min(pos(:,1));
[right,indi]=max(pos(:,1));
right=right+pos(indi,3);
bot=pos(1,2)+pos(1,4);
wid=right-left;
ta=axes('Position',[left bot wid 0.1]);
th=text(0.5,0.2,tstr,'FontSize',fs,'HorizontalAlignment','center');
set(ta,'visible','off')

varns={th,ta};
varargout=varns(1:nargout);
