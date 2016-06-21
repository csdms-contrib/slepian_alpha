function [ax,xl,yl]=xtraxis(ah,xti,xtil,xlab,yti,ytil,ylab)
% [ax,xl,yl]=xtraxis(aha,xti,xtil,xlab,yti,ytil,ylab)
%
% Creates an extra axis on top of the current one
% 
% INPUT
%
% ah        Axis handle (scalar) (default: gca)
% xti       Major positions to be labeled on the new x-axis
% xtil      Labels for xti to be put on new x-axis
% xlab      Label string for the new x-axis itself 
% yti       Positions to be labeled on both old and new y-axis
% ytil      Labels for yti to be put on new y-axis
% ylab      Label string for the new y-axis itself
%
% OUTPUT
%
% ax        Handle to the extra axis
% xl,yl     Handles to the x and y labels 
%
% Best usage for the logarithmic case is shown in, e.g., MERMAID04,
% COHANAL2, EGGERS2, etc. where it's all different 
% xtraxis(ah,lambda,lambda) is key, and above all, need 'xdir' 'rev' and
% specify 'xlim' (using the parent values) again, at the end.
% A recent example of the linear case is in SWREGIONS2D. 
%
% See LAXIS, XTRAXIS1D, XTRAXIS2D
%
% Last modified by fjsimons-at-alum.mit.edu, 06/08/2015

defval('ah',gca)
defval('xti',[])
defval('xtil',[])
defval('yti',[])
defval('ytil',[])
defval('xlab',[])
defval('ylab',[])
defval('xl',[])
defval('yl',[])

% Make sure ticks are unique
[xti,j]=unique(xti);
if ~isstr(xtil) && ~isempty(xtil)
  xtil=xtil(j);
end

[yti,j]=unique(yti);
if ~isstr(ytil) & ~isempty(ytil)
  ytil=ytil(j);
end

% Create the axis
[ax,axl,loc]=laxis(ah,0,0);

% Put on the same limits as the old one
xll=xlim(ah);
yll=ylim(ah);
set(ah,'box','off')
set(ax,'XTick',xti,'XTickLabel',xtil,...
       'xlim',xll,'ylim',yll,...
       'ytick',yti,'ytickLabel',ytil,...
       'Xgrid','off','Ygrid','off',...
       'XaxisLabel','top','YaxisLabel','right',...
       'XScale',get(ah,'XScale'),'YScale',get(ah,'YScale'),...
       'box','off','Visible','on','Color','none');
% For some reason these are tricky when there's auto
set(ax,'DataAspectRatio',get(ah,'DataAspectRatio'),...
       'PlotBoxAspectRatio',get(ah,'PlotBoxAspectRatio'),...
       'DataAspectRatioMode',get(ah,'DataAspectRatioMode'),...
       'PlotBoxAspectRatioMode',get(ah,'PlotBoxAspectRatioMode'))
axes(ax)
xl=xlabel(xlab);
yl=ylabel(ylab);
