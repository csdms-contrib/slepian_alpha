function cbarticks(cb,cax,tint,pos)
% CBARTICKS(cb,cax,tint,pos)
%
% Transforms axis tick marks and labels to other values, e.g. for a
% colorbar that indexes an IMAGE type figure, so that the ticks, which
% would otherwise be marking color indices, agree with the actual values
% of what's being plotted in the figure and for which the axis is a
% colorbar. 
%
% INPUT:
%
% cb       Handle to an axis, especially, a colorbar from ADDCB
% cax      Minimum and maximum of the color range
% tint     If a single number: tick INTERVAL, cax(1):tint:cax(2)
%          If two or more numbers: actual ticks themselves
% pos      'hor' or 'vert' position of the axis handle cb
%
% SEE ALSO:
%
% ADDCB, IMAGEF, IMAGEFDIR, IMAGEFNAN
% 
% Last modified by charig-at-princeton.edu, 02/05/2016
% Last modified by fjsimons-at-alum.mit.edu, 07/26/2017

defval('tint',20)
defval('pos','hor')

if length(tint)==1
  % Make sure the last one makes it regardless of tint
  tix=unique([cax(1):tint:cax(2) cax(2)]);
else
  tix=tint;
end

if strcmp(pos,'hor')
  set(cb,'Xtick',scale(tix,get(cb,'xlim')))
  set(cb,'XtickLabel',tix)
elseif strcmp(pos,'vert')
  set(cb,'Ytick',scale(tix,get(cb,'ylim')))
  set(cb,'YtickLabel',tix)
else
  error('Specify valid orientation')
end


