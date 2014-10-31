function [ax,xl,yl]=xtraxis(ah,xti,xtil,xlab,yti,ytil,ylab)
% [ax,xl,yl]=xtraxis(aha,xti,xtil,xlab,yti,ytil,ylab)
%
% Creates an extra axis on top of the current one
% 
% INPUT
%
% ah        Axis handle (scalar) (default: gca)
% xti       Positions to be labeled on both old and new x-axis
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
% Best usage for the logarithmic case is shown in MERMAID04 and COHANAL
% where it's all different and a recent example of the linear case is in
% SWREGIONS2d. Above all, need 'xdir' 'rev' and specify 'xlim' again.
%
% See LAXIS, XTRAXIS1D, XTRAXIS2D
%
% Last modified by fjsimons-at-alum.mit.edu, 04/01/2010

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
set(ax,'xtick',xti,'xtickl',xtil,'xlim',xll,'ylim',yll,...
       'ytick',yti,'ytickl',ytil,...
       'Visible','on','Color','none','Xgrid','off','Ygrid','off',...
       'XaxisL','top','DataAspectRatio',get(ah,'DataAspectRatio'),...
       'PlotBoxAspectRatio',get(ah,'PlotBoxAspectRatio'),...
       'PlotBoxAspectRatioMode',get(ah,'PlotBoxAspectRatioMode'),...
       'DataAspectRatioMode',get(ah,'DataAspectRatioMode'),...
       'XScale',get(ah,'XScale'),...
       'YScale',get(ah,'YScale'),...
       'YaxisL','right','box','off');
axes(ax)
xl=xlabel(xlab);
yl=ylabel(ylab);

