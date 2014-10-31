function varargout=supertit(ah,tstr,fs)
% [H,ta]=SUPERTIT(ah,tstr,fs)
%
% Takes axis ah and puts a title tstr spanning over the panels in
% question. These ah need to refer to plots on the same window row,
% and need to have the same height. 
%
% INPUT:
%
% ah        Axis handles
% tstr      Title string
% fs        Fontsize [default 15]
%
% OUTPUT:
%
% H         Handle to the supertitle
% ta        Handle to the invisible axis containing it
%
%
% EXAMPLE:
%
% a=krijetem([231 232 233]);
% b=supertit(a([1 2]),'This is the supertitle');
%
% Last modified by fjsimons-at-alum.mit.edu, 08/15/2007

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
th=text(0.5,0.2,tstr,'FontS',fs,'HorizontalAlignment','center');
set(ta,'visible','off')

varns={th,ta};
varargout=varns(1:nargout);
