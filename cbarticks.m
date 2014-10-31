function cbarticks(cb,cax,tint,pos)
% CBARTICKS(cb,cax,tint,pos)
%
% Transforms the tick marks and labels to agree with
% what's plotted in the figure.
%
% INPUT:
%
% cb       Handle to the colorbar
% cax      Minimum and maximum of the color range
% tint     Tick interval
% pos      'hor' or 'vert'
%
% SEE ALSO:
%
% ADDCB, IMAGEF, IMAGEFDIR, IMAGEFNAN
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/19/2008

defval('tint',20)
defval('pos','hor')

% Make sure the last one makes it regardless of tint
tix=unique([cax(1):tint:cax(2) cax(2)]);

if strcmp(pos,'hor')
  set(cb,'Xtick',scale(tix,get(cb,'xlim')))
  set(cb,'XtickL',tix)
elseif strcmp(pos,'vert')
  set(cb,'Ytick',scale(tix,get(cb,'ylim')))
  set(cb,'YtickL',tix)
else
  error('Specify valid orientation')
end


